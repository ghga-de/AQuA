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
import re
from collections import Counter

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


def _first_value(record: dict, keys: tuple[str, ...]) -> str:
    """Return the first non-empty scalar/list value from a metadata record."""
    for key in keys:
        value = record.get(key)
        if value is None or value == "":
            continue
        values = _as_list(value)
        return ";".join(values) if values else str(value).strip()
    return ""


def get_analysis_tool(analysis: dict, analysis_method_record: dict) -> str:
    """Extract a human-readable analysis tool or workflow name when available."""
    return _first_value(
        analysis,
        (
            "analysis_tool",
            "tool",
            "tools",
            "software",
            "software_name",
            "workflow",
            "workflow_name",
            "pipeline",
            "pipeline_name",
        ),
    ) or _first_value(
        analysis_method_record,
        (
            "analysis_tool",
            "tool",
            "tools",
            "software",
            "software_name",
            "workflow",
            "workflow_name",
            "pipeline",
            "pipeline_name",
            "name",
        ),
    )


def get_analysis_genome(analysis: dict, analysis_method_record: dict) -> str:
    """Extract reference genome / assembly metadata when the submission provides it."""
    return _first_value(
        analysis,
        (
            "analysis_genome",
            "reference_genome",
            "reference_assembly",
            "genome_assembly",
            "genome",
            "assembly",
            "reference",
        ),
    ) or _first_value(
        analysis_method_record,
        (
            "analysis_genome",
            "reference_genome",
            "reference_assembly",
            "genome_assembly",
            "genome",
            "assembly",
            "reference",
        ),
    )


def get_target_bed(analysis: dict, analysis_method_record: dict) -> str:
    """Extract target BED / interval-set metadata when available."""
    return _first_value(
        analysis,
        (
            "target_bed",
            "bed",
            "bed_file",
            "capture_bed",
            "target_intervals",
            "target_interval_file",
            "intervals",
            "panel_bed",
        ),
    ) or _first_value(
        analysis_method_record,
        (
            "target_bed",
            "bed",
            "bed_file",
            "capture_bed",
            "target_intervals",
            "target_interval_file",
            "intervals",
            "panel_bed",
        ),
    )


def get_alignment_status(analysis: dict, file_record: dict) -> str:
    """Return aligned/unaligned status if metadata names it explicitly."""
    raw = _first_value(
        file_record,
        ("alignment_status", "read_alignment_status", "aligned", "is_aligned"),
    ) or _first_value(
        analysis,
        ("alignment_status", "read_alignment_status", "aligned", "is_aligned"),
    )
    value = raw.strip().lower()
    if value in {"true", "yes", "aligned"}:
        return "aligned"
    if value in {"false", "no", "unaligned", "unmapped"}:
        return "unaligned"
    return value


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
    Returns a dict with keys: fastq_1, fastq_2, bam, bai, cram, crai, vcf, bed, other.
    """
    buckets: dict[str, list] = {
        k: []
        for k in (
            "fastq_1",
            "fastq_2",
            "bam",
            "bai",
            "cram",
            "crai",
            "vcf",
            "bed",
            "other",
        )
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
        elif fmt in {"BED", "BED12"} or low.endswith((".bed", ".bed.gz")):
            buckets["bed"].append(path)
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

    file_rows.extend(
        pair_fastqs(buckets["fastq_1"], buckets["fastq_2"], single_end, warnings, context)
    )
    target_bed = buckets["bed"][0] if buckets["bed"] else ""
    if len(buckets["bed"]) > 1:
        warnings.append(
            f"[{context}] multiple BED files found; using {target_bed!r} as target_bed "
            f"and ignoring {buckets['bed'][1:]}."
        )

    def _with_target_bed(row: dict) -> dict:
        if target_bed:
            return {**row, "target_bed": target_bed}
        return row

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
        file_rows.append(
            _with_target_bed({"bam": bam, "bai": bai, "single_end": single_end})
        )

    for cram in buckets["cram"]:
        crai = crai_by_stem.get(_stem(cram), "")
        file_rows.append(
            _with_target_bed({"cram": cram, "crai": crai, "single_end": single_end})
        )

    for vcf in buckets["vcf"]:
        file_rows.append(_with_target_bed({"vcf": vcf, "single_end": single_end}))

    if buckets["other"] and not file_rows:
        file_rows.append(
            {
                "data_files": ";".join(buckets["other"]),
                "target_bed": target_bed,
                "single_end": single_end,
            }
        )

    if target_bed and not file_rows:
        file_rows.append(
            {"data_files": target_bed, "target_bed": target_bed, "single_end": single_end}
        )

    # Return None when no files were found so the caller can skip the row
    # rather than writing a metadata-only row with no usable file inputs.
    if not file_rows:
        return None

    return file_rows


def pair_fastqs(
    fastq_1: list[str],
    fastq_2: list[str],
    single_end: str,
    warnings: list[str],
    context: str = "",
) -> list[dict]:
    """Pair FASTQs by normalized mate basename instead of only sorted position."""

    def _pair_key(path: str) -> str:
        base = os.path.basename(path).lower()
        base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base)
        patterns = (
            r"([._-])r[12]([._-]?\d+)?$",
            r"([._-])[12]$",
            r"([._-])r[12]([._-])",
            r"([._-])[12]([._-])",
        )
        for pattern in patterns:
            stripped = re.sub(pattern, r"\1", base)
            if stripped != base:
                return stripped.rstrip("._-")
        return base

    r1_by_key: dict[str, list[str]] = {}
    r2_by_key: dict[str, list[str]] = {}
    for path in fastq_1:
        r1_by_key.setdefault(_pair_key(path), []).append(path)
    for path in fastq_2:
        r2_by_key.setdefault(_pair_key(path), []).append(path)

    rows = []
    all_keys = sorted(set(r1_by_key) | set(r2_by_key))
    for key in all_keys:
        r1s = sorted(r1_by_key.get(key, []))
        r2s = sorted(r2_by_key.get(key, []))
        if len(r1s) > 1 or len(r2s) > 1:
            warnings.append(
                f"[{context}] ambiguous FASTQ mate group {key!r}: "
                f"R1={r1s}, R2={r2s}; pairing by sorted order within the group."
            )
        max_len = max(len(r1s), len(r2s))
        for idx in range(max_len):
            f1 = r1s[idx] if idx < len(r1s) else ""
            f2 = r2s[idx] if idx < len(r2s) else ""
            if not f1:
                warnings.append(
                    f"[{context}] FASTQ mate group {key!r} has fastq_2={f2!r} "
                    "without a matching fastq_1."
                )
            if not f2 and single_end == "false":
                warnings.append(
                    f"[{context}] PE experiment is missing fastq_2 for fastq_1={f1!r}; "
                    "single_end kept as 'false' from metadata — check for missing R2 file."
                )
            rows.append(
                {"fastq_1": f1 or "", "fastq_2": f2 or "", "single_end": single_end}
            )

    return rows


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

        ### Collect files in provenance groups instead of flattening the sample.
        ### Raw files are grouped by experiment; processed files are grouped by
        ### analysis so BAM/CRAM/VCF rows keep their own analysis metadata.
        file_groups: list[dict] = []
        seen_raw_names: set[str] = set()
        seen_process_names: set[str] = set()
        seen_analysis_aliases: set[str] = set()

        for experiment in idx["sample_to_experiments"].get(sample_alias, []):
            exp_alias = experiment.get("alias", "")
            em_alias = experiment.get("experiment_method", "")
            exp_method = idx["exp_methods"].get(em_alias, {})
            exp_method_name = get_method(exp_method)  # captured for this experiment
            single_end = get_single_end(exp_method)  # captured for this experiment

            ### Raw files – tagged with this experiment's provenance.
            raw_files = []
            for rdf in idx["exp_to_rdfs"].get(exp_alias, []):
                name = rdf.get("name") or rdf.get("alias") or ""
                if name and name not in seen_raw_names:
                    seen_raw_names.add(name)
                    raw_files.append(rdf)

            raw_analysis_methods = []
            for analysis in idx["exp_to_analyses"].get(exp_alias, []):
                am_record = idx["analysis_methods"].get(
                    analysis.get("analysis_method", ""), {}
                )
                raw_analysis_methods.append(get_analysis_method(am_record))

            if raw_files:
                raw_counts = Counter(raw_analysis_methods)
                raw_analysis_method = (
                    raw_counts.most_common(1)[0][0] if raw_counts else "unknown"
                )
                file_groups.append(
                    {
                        "files": raw_files,
                        "experiment_method": exp_method_name,
                        "analysis_method": raw_analysis_method,
                        "analysis_tool": "",
                        "analysis_genome": "",
                        "target_bed": "",
                        "alignment_status": "",
                        "single_end": single_end,
                        "context": f"{sample_alias}/{exp_alias}/raw",
                    }
                )

            ### Process files – tagged with their specific analysis provenance
            for analysis in idx["exp_to_analyses"].get(exp_alias, []):
                analysis_alias = analysis.get("alias", "")
                if analysis_alias in seen_analysis_aliases:
                    continue
                seen_analysis_aliases.add(analysis_alias)

                am_alias = analysis.get("analysis_method", "")
                am_record = idx["analysis_methods"].get(am_alias, {})
                ana_method_name = get_analysis_method(
                    am_record
                )  # captured for this analysis
                process_files = []
                for pdf in idx["analysis_to_pdfs"].get(analysis_alias, []):
                    name = pdf.get("name") or pdf.get("alias") or ""
                    if name and name not in seen_process_names:
                        seen_process_names.add(name)
                        process_files.append(pdf)
                if not process_files:
                    continue

                alignment_statuses = [
                    get_alignment_status(analysis, pdf) for pdf in process_files
                ]
                alignment_statuses = [status for status in alignment_statuses if status]
                alignment_status = (
                    Counter(alignment_statuses).most_common(1)[0][0]
                    if alignment_statuses
                    else ""
                )
                file_groups.append(
                    {
                        "files": process_files,
                        "experiment_method": exp_method_name,
                        "analysis_method": ana_method_name,
                        "analysis_tool": get_analysis_tool(analysis, am_record),
                        "analysis_genome": get_analysis_genome(analysis, am_record),
                        "target_bed": get_target_bed(analysis, am_record),
                        "alignment_status": alignment_status,
                        "single_end": single_end,
                        "context": f"{sample_alias}/{analysis_alias}",
                    }
                )

        if not file_groups:
            warnings.append(
                f"[{sample_alias}/no-files] no files found in any bucket — "
                "skipping row to avoid writing a metadata-only entry."
            )
            continue

        sample_base = {
            "sample": sample_alias,
            "individual_id": individual_id,
            "sex": sex,
            "status": status,
            "phenotype": phenotype,
            "sample_type": sample.get("type", ""),
            "disease_status": sample.get("disease_or_healthy", ""),
            "case_control_status": sample.get("case_control_status", ""),
            "tissue": sample.get("biospecimen_tissue_term", ""),
        }

        lane_num = 1
        for group in file_groups:
            buckets = classify_files(group["files"], input_directory)
            file_rows = buckets_to_rows(
                buckets, group["single_end"], warnings, group["context"]
            )

            if file_rows is None:
                warnings.append(
                    f"[{group['context']}] no files found in any bucket — "
                    "skipping row to avoid writing a metadata-only entry."
                )
                continue

            base = {
                **sample_base,
                "experiment_method": group["experiment_method"],
                "analysis_method": group["analysis_method"],
                "analysis_tool": group["analysis_tool"],
                "analysis_genome": group["analysis_genome"],
                "target_bed": group["target_bed"],
                "alignment_status": group["alignment_status"],
            }

            for frow in file_rows:
                row = {**base, "lane": f"L{lane_num:03d}", **frow}
                rows.append(row)
                lane_num += 1

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
    "analysis_tool",
    "analysis_genome",
    "target_bed",
    "alignment_status",
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
