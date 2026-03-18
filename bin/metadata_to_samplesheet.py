#!/usr/bin/env python3
"""
Convert a GHGA metadata JSON (submission format) to a pipeline samplesheet CSV.

GHGA metadata model reference: https://docs.ghga.de/metadata/

Key entities and how they map to samplesheet columns:

  Individual        --> individual_id, sex, phenotype
  Sample            --> sample, status, case_control_status, tissue, sample_type, disease_status
  ExperimentMethod  --> experiment_method (normalised: "wgs", "pacbio", …),
                        experiment_method_alias (original alias, for traceability),
                        single_end
  Experiment        --> links Sample ↔ ExperimentMethod ↔ ResearchDataFiles
  ResearchDataFile  --> fastq_1, fastq_2  (raw FASTQ/FAST5; sorted by technical_replicate)
  Analysis          --> analysis_method (normalised), analysis_method_alias;
                        bridges Experiments ↔ ProcessDataFiles
  ProcessDataFile   --> bam/bai, cram/crai, vcf  (linked via Analysis)

ProcessDataFile linkage chain:
  ProcessDataFile.analysis --> Analysis.alias
  Analysis.research_data_files --> ResearchDataFile.alias (list)
  ResearchDataFile.experiments --> Experiment.alias (list)
  Experiment.sample --> Sample.alias

Usage
-----
  python metadata_to_samplesheet.py metadata.json
  python metadata_to_samplesheet.py metadata.json --output samplesheet.csv
  python metadata_to_samplesheet.py metadata.json --input_directory /data/files
"""

import argparse
import csv
import json
import os
from collections import Counter
from itertools import zip_longest
from pathlib import Path

##########################
### Method resolution  ###
##########################

LIBRARY_TYPE_MAP = {
    "WGS": "wgs",
    "WXS": "wes",
    "WCS": "wes",
    "TOTAL_RNA": "rna",
    "M_RNA": "rna",
    "MI_RNA": "smrna",
    "NC_RNA": "rna",
    "ATAC": "atac",
    "METHYLATION": "methylseq",
    "CHIP_SEQ": "chip",
    "OTHER": "wgs",  # fallback for custom panels (e.g. cfDNA)
}

NANOPORE_INSTRUMENTS = {"MINION", "GRIDION", "PROMETHION"}
PACBIO_INSTRUMENTS = {"PACBIO_RS", "PACBIO_RS_II", "SEQUEL", "SEQUEL_II", "SEQUEL_IIE"}


def get_method(exp_method: dict) -> str:
    """Resolve pipeline method from instrument_model (takes priority) or library_type."""
    instrument = exp_method.get("instrument_model", "").upper()
    if instrument in NANOPORE_INSTRUMENTS:
        return "nanopore"
    if instrument in PACBIO_INSTRUMENTS:
        return "pacbio"
    library_type = exp_method.get("library_type", "").strip().upper()
    return LIBRARY_TYPE_MAP.get(library_type, "unknown")


def get_analysis_method(analysis_method_record: dict) -> str:
    """
    Resolve a normalised pipeline method name from an AnalysisMethod record.
    Uses the 'type' field, mapped through LIBRARY_TYPE_MAP for consistency with
    get_method (e.g. 'cfDNA' stays as-is when not in the map, lowercased).
    Falls back to 'unknown' when the record is absent.
    """
    if not analysis_method_record:
        return "unknown"
    ana_type = analysis_method_record.get("type", "")
    return LIBRARY_TYPE_MAP.get(ana_type, ana_type.lower() or "unknown")


def get_single_end(exp_method: dict) -> str:
    """
    Derive single_end flag from ExperimentMethod.sequencing_layout.
    'SE' --> 'true'; 'PE' (or anything else) --> 'false'.
    """
    layout = exp_method.get("sequencing_layout", "").upper()
    return "true" if layout == "SE" else "false"


##########################
###   Index Builders   ###
##########################


def _as_list(value) -> list:
    """Normalise a field that may be a list, a comma-separated string, or a scalar."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v).strip() for v in value if v is not None]
    return [v.strip() for v in str(value).split(",") if v.strip()]


def build_indices(metadata: dict) -> dict:
    """
    Pre-build all lookup dictionaries needed for samplesheet construction.
    Returns a dict of named indices.
    """
    ### individuals: alias --> record
    individuals = {
        ind["alias"]: ind for ind in metadata.get("individuals", []) if ind.get("alias")
    }

    ### experiment_methods: alias --> record
    exp_methods = {
        em["alias"]: em
        for em in metadata.get("experiment_methods", [])
        if em.get("alias")
    }

    ### experiments: alias --> record
    experiments = {
        exp["alias"]: exp for exp in metadata.get("experiments", []) if exp.get("alias")
    }

    ### sample --> experiments (all experiments; a sample may have multiple sequencing runs)
    sample_to_experiments: dict[str, list] = {}
    for exp in metadata.get("experiments", []):
        sample_ref = exp.get("sample")
        if sample_ref:
            sample_to_experiments.setdefault(sample_ref, []).append(exp)

    ### experiment --> ResearchDataFiles (list)
    exp_to_rdfs: dict[str, list] = {}
    for rdf in metadata.get("research_data_files", []):
        for exp_alias in _as_list(rdf.get("experiments")):
            exp_to_rdfs.setdefault(exp_alias, []).append(rdf)

    ### ResearchDataFile alias --> experiment aliases
    rdf_alias_to_experiments: dict[str, list] = {}
    for rdf in metadata.get("research_data_files", []):
        rdf_alias = rdf.get("alias")
        if rdf_alias:
            rdf_alias_to_experiments[rdf_alias] = _as_list(rdf.get("experiments"))

    ### analysis_methods: alias --> record
    analysis_methods = {
        am["alias"]: am
        for am in metadata.get("analysis_methods", [])
        if am.get("alias")
    }

    ### analyses: alias --> record
    analyses = {
        ana["alias"]: ana for ana in metadata.get("analyses", []) if ana.get("alias")
    }

    ### experiment --> analyses (all analyses linked via RDF membership;
    ### an experiment may be re-analysed multiple times)
    ### Analysis.research_data_files --> list of RDF aliases --> each RDF --> experiments
    exp_to_analyses: dict[str, list] = {}
    for ana in metadata.get("analyses", []):
        for rdf_alias in _as_list(ana.get("research_data_files")):
            for exp_alias in rdf_alias_to_experiments.get(rdf_alias, []):
                if ana not in exp_to_analyses.get(exp_alias, []):
                    exp_to_analyses.setdefault(exp_alias, []).append(ana)

    ### ProcessDataFiles: analysis alias --> list of process files
    ### (ProcessDataFile links to Analysis, not directly to Experiment)
    analysis_to_pdfs: dict[str, list] = {}
    for pdf in metadata.get("process_data_files", []):
        ana_ref = pdf.get("analysis")
        if ana_ref:
            analysis_to_pdfs.setdefault(ana_ref, []).append(pdf)

    return {
        "individuals": individuals,
        "exp_methods": exp_methods,
        "experiments": experiments,
        "sample_to_experiments": sample_to_experiments,
        "exp_to_rdfs": exp_to_rdfs,
        "exp_to_analyses": exp_to_analyses,
        "analysis_to_pdfs": analysis_to_pdfs,
        "analyses": analyses,
        "analysis_methods": analysis_methods,
    }


#############################
#### File classification ####
#############################


def classify_files(files: list[dict], input_directory: str) -> dict[str, list]:
    """
    Sort a flat list of file records into typed buckets.
    Falls back to filename extension when format field is missing/ambiguous.
    Returns a dict with keys: fastq_1, fastq_2, bam, bai, cram, crai, vcf, other.
    """
    buckets: dict[str, list] = {
        k: []
        for k in ("fastq_1", "fastq_2", "bam", "bai", "cram", "crai", "vcf", "other")
    }
    seen: set = set()

    for f in files:
        name = f.get("name") or f.get("alias") or ""
        if not name or name in seen:
            continue
        seen.add(name)

        path = os.path.join(input_directory, name) if input_directory else name
        fmt = str(f.get("format", "")).upper()
        low = name.lower()

        if fmt == "FASTQ" or low.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
            # Only schema-allowed FASTQ format goes into fastq buckets.
            # Use technical_replicate when it explicitly says 2 (definitive R2).
            # Otherwise fall back to filename patterns: _R2, _2.fastq, etc.
            # This handles cfDNA panels where replicate tracks library (both=1),
            # as well as standard WGS/WES where R1=1, R2=2.
            replicate = f.get("technical_replicate")
            is_r2_by_name = "_r2" in low or "_2.f" in low or ".r2." in low
            if replicate == 2 or is_r2_by_name:
                buckets["fastq_2"].append(path)
            else:
                buckets["fastq_1"].append(path)

        elif fmt == "BAI" or low.endswith(".bai"):
            buckets["bai"].append(path)
        elif fmt == "BAM" or low.endswith(".bam"):
            # SAM intentionally excluded: not a schema-allowed column
            buckets["bam"].append(path)
        elif fmt == "CRAI" or low.endswith(".crai"):
            buckets["crai"].append(path)
        elif fmt == "CRAM" or low.endswith(".cram"):
            buckets["cram"].append(path)
        elif fmt == "VCF" or low.endswith((".vcf", ".vcf.gz")):
            # BCF intentionally excluded: not a schema-allowed column
            buckets["vcf"].append(path)
        else:
            # FAST5, FASTA, UBAM, SAM, BCF and any other non-schema formats
            buckets["other"].append(path)

    for key in buckets:
        buckets[key].sort()

    return buckets


def buckets_to_rows(
    buckets: dict[str, list],
    single_end: str,
    warnings: list[str],
    context: str = "",
) -> list[dict] | None:
    """
    Expand classified file buckets into one or more file-row dicts.
    Each row carries exactly one file type (FASTQ pair, BAM+BAI, CRAM+CRAI, or VCF).

    Returns None when no files exist in any bucket
    Single_end is derived from metadata (sequencing_layout).
    When a paired-end experiment is missing its R2 file a warning is appended to `warnings`
    """
    file_rows = []

    for f1, f2 in zip_longest(buckets["fastq_1"], buckets["fastq_2"]):
        if not f2 and single_end == "false":
            warnings.append(
                f"[{context}] PE experiment is missing fastq_2 for fastq_1={f1!r}; "
                "single_end kept as 'false' from metadata — check for missing R2 file."
            )
        file_rows.append(
            {"fastq_1": f1 or "", "fastq_2": f2 or "", "single_end": single_end}
        )

    def _stem(path: str) -> str:
        """Strip directory and all extensions to get a bare basename for index matching."""
        base = os.path.basename(path)
        # Strip known multi-part extensions first (.vcf.gz, .fastq.gz, etc.)
        for ext in (".bam", ".cram", ".bai", ".crai"):
            if base.endswith(ext):
                return base[: -len(ext)]
        # Fall back to stripping whatever os.path.splitext finds
        return os.path.splitext(base)[0]

    # Build stem --> path maps for index files so BAM/CRAM can look up their index
    bai_by_stem = {_stem(p): p for p in buckets["bai"]}
    crai_by_stem = {_stem(p): p for p in buckets["crai"]}

    for bam in buckets["bam"]:
        bai = bai_by_stem.get(_stem(bam), "")
        file_rows.append({"bam": bam, "bai": bai, "single_end": single_end})

    for cram in buckets["cram"]:
        crai = crai_by_stem.get(_stem(cram), "")
        file_rows.append({"cram": cram, "crai": crai, "single_end": single_end})

    for vcf in buckets["vcf"]:
        file_rows.append({"vcf": vcf, "single_end": single_end})

    if buckets["other"] and not file_rows:
        file_rows.append(
            {"data_files": ";".join(buckets["other"]), "single_end": single_end}
        )

    # Return None when no files were found so the caller can skip the row
    # rather than writing a metadata-only row with no usable file inputs.
    if not file_rows:
        return None

    return file_rows


#########################
#####  Main parser  #####
#########################


def parse_metadata(metadata_path: str, input_directory: str) -> list[dict]:
    with open(metadata_path) as f:
        metadata = json.load(f)

    idx = build_indices(metadata)
    warnings: list[str] = []

    rows = []
    for sample in metadata.get("samples", []):
        sample_alias = sample.get("alias", "")
        if not sample_alias:
            continue

        ### Individual
        ind_ref = sample.get("individual", "")
        individual = idx["individuals"].get(ind_ref, {})
        individual_id = individual.get("alias") or ind_ref or ""

        sex_raw = individual.get("sex", "")
        sex = str(sex_raw).strip() if sex_raw else "NA"

        phenotype_raw = individual.get("phenotypic_features_terms", [])
        phenotype = ";".join(_as_list(phenotype_raw)) if phenotype_raw else ""

        ### Sample metadata
        disease_raw = str(sample.get("disease_or_healthy", "")).lower()
        status = 1 if ("tumor" in disease_raw or "disease" in disease_raw) else 0

        ### Collect all files for this sample across all experiments and analyses.
        ### Each file is tagged with its own provenance (exp_method_name,
        ### ana_method_name, single_end) captured at the point of collection so
        ### earlier files are never overwritten by later loop iterations.
        ### Files are deduplicated by name across the whole sample.

        # Each entry: (file_record, exp_method_name, ana_method_name, single_end)
        sample_files: list[tuple[dict, str, str, str]] = []
        seen_names: set[str] = set()

        for experiment in idx["sample_to_experiments"].get(sample_alias, []):
            exp_alias = experiment.get("alias", "")
            em_alias = experiment.get("experiment_method", "")
            exp_method = idx["exp_methods"].get(em_alias, {})
            exp_method_name = get_method(exp_method)  # captured for this experiment
            single_end = get_single_end(exp_method)  # captured for this experiment

            ### Raw files – tagged with this experiment's provenance
            for rdf in idx["exp_to_rdfs"].get(exp_alias, []):
                name = rdf.get("name") or rdf.get("alias") or ""
                if name and name not in seen_names:
                    seen_names.add(name)
                    # ana_method_name resolved below once all analyses are known
                    sample_files.append((rdf, exp_method_name, None, single_end))

            ### Process files – tagged with their specific analysis provenance
            for analysis in idx["exp_to_analyses"].get(exp_alias, []):
                analysis_alias = analysis.get("alias", "")
                am_alias = analysis.get("analysis_method", "")
                am_record = idx["analysis_methods"].get(am_alias, {})
                ana_method_name = get_analysis_method(
                    am_record
                )  # captured for this analysis

                for pdf in idx["analysis_to_pdfs"].get(analysis_alias, []):
                    name = pdf.get("name") or pdf.get("alias") or ""
                    if name and name not in seen_names:
                        seen_names.add(name)
                        sample_files.append(
                            (pdf, exp_method_name, ana_method_name, single_end)
                        )

        if not sample_files:
            warnings.append(
                f"[{sample_alias}/no-files] no files found in any bucket — "
                "skipping row to avoid writing a metadata-only entry."
            )
            continue

        ### Resolve ana_method_name for raw files (am=None) using the most
        ### common analysis method seen for this sample; ties broken by first seen.
        ana_counts = Counter(am for _, _, am, _ in sample_files if am is not None)
        dominant_ana_method = (
            ana_counts.most_common(1)[0][0] if ana_counts else "unknown"
        )

        resolved_files = [
            (f, em, am if am is not None else dominant_ana_method, se)
            for f, em, am, se in sample_files
        ]

        ### Use provenance from the first file as the row-level method values.
        ### Within a well-formed submission all files on a sample share the same
        ### experiment/analysis method; if they differ the first value is used
        ### and the per-file tags in resolved_files remain available for future use.
        _, row_exp_method, row_ana_method, row_single_end = resolved_files[0]

        ### Classify and emit once — single lane counter for the whole sample
        all_file_records = [f for f, _, _, _ in resolved_files]
        context = sample_alias
        buckets = classify_files(all_file_records, input_directory)
        file_rows = buckets_to_rows(buckets, row_single_end, warnings, context)

        if file_rows is None:
            warnings.append(
                f"[{context}] no files found in any bucket — "
                "skipping row to avoid writing a metadata-only entry."
            )
            continue

        base = {
            "sample": sample_alias,
            "individual_id": individual_id,
            "sex": sex,
            "status": status,
            "phenotype": phenotype,
            "sample_type": sample.get("type", ""),
            "disease_status": sample.get("disease_or_healthy", ""),
            "case_control_status": sample.get("case_control_status", ""),
            "tissue": sample.get("biospecimen_tissue_term", ""),
            "experiment_method": row_exp_method,
            "analysis_method": row_ana_method,
        }

        for lane_num, frow in enumerate(file_rows, start=1):
            row = {**base, "lane": f"L{lane_num:03d}", **frow}
            rows.append(row)

    if warnings:
        import sys

        print(f"\n{len(warnings)} validation warning(s):", file=sys.stderr)
        for w in warnings:
            print(f"  WARNING: {w}", file=sys.stderr)

    return rows


####################
#####  Writer  #####
####################

ALL_COLUMNS = [
    "sample",
    "lane",
    "individual_id",
    "sex",
    "status",
    "phenotype",
    "sample_type",
    "disease_status",
    "case_control_status",
    "tissue",
    "experiment_method",
    "analysis_method",
    "fastq_1",
    "fastq_2",
    "single_end",
    "bam",
    "bai",
    "cram",
    "crai",
    "vcf",
    "data_files",
]

# Columns that must be written as lowercase "true"/"false" strings to match
# checked-in samplesheets and satisfy the schema's "type": "boolean" constraint.
_BOOL_COLUMNS = {"single_end"}


def write_samplesheet(rows: list[dict], output_path: str):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ALL_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = {}
            for col in ALL_COLUMNS:
                val = row.get(col, "")
                if col in _BOOL_COLUMNS and val != "":
                    val = "true" if str(val).lower() == "true" else "false"
                out[col] = val
            writer.writerow(out)


#############
#### CLI ####
#############


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert GHGA metadata JSON to pipeline samplesheet CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("metadata", help="Path to GHGA metadata.json")
    parser.add_argument(
        "--output",
        default="samplesheet.csv",
        help="Output CSV path (default: samplesheet.csv)",
    )
    parser.add_argument(
        "--input_directory",
        default="",
        help="Optional prefix directory prepended to all file names in the output",
    )
    args = parser.parse_args()

    rows = parse_metadata(args.metadata, args.input_directory)
    write_samplesheet(rows, args.output)
    print(f"Written {len(rows)} rows to {args.output}")
