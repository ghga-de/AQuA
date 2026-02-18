#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path
import pandas as pd
from itertools import zip_longest

def sanitize(s: str) -> str:
    if not s:
        return ""
    return str(s).strip().replace(" ", "_")

def main(input_json: Path, input_directory: str, output_csv: str):
    with open(input_json, "r") as json_file:
        data = json.load(json_file)

    individuals = {}
    for ind in data.get("individuals", []):
        individuals[ind.get("alias")] = ind
        if ind.get("accession"):
            individuals[ind.get("accession")] = ind

    analyses = {}
    for ana in data.get("analyses", []):
        analyses[ana.get("alias")] = ana
        if ana.get("accession"):
            analyses[ana.get("accession")] = ana

    data_files = {}
    for k in ["research_data_files", "process_data_files", "processed_data_files"]:
        for f in data.get(k, []):
            data_files[f.get("alias")] = f
            if f.get("accession"):
                data_files[f.get("accession")] = f

    exp_to_sample = {}
    sample_to_exp_method = {}
    for exp in data.get("experiments", []):
        exp_to_sample[exp.get("alias")] = exp.get("sample")
        if exp.get("accession"):
            exp_to_sample[exp.get("accession")] = exp.get("sample")
            
        sample_ref = exp.get("sample")
        if sample_ref:
            sample_to_exp_method[sample_ref] = exp.get("experiment_method", "")

    sample_to_files = {}
    for k in ["research_data_files", "process_data_files", "processed_data_files"]:
        for f in data.get(k, []):
            exp_refs = f.get("experiments", [])
            if isinstance(exp_refs, str):
                exp_refs = [e.strip() for e in exp_refs.split(",") if e.strip()]
            elif not isinstance(exp_refs, list):
                exp_refs = [exp_refs]
                
            for exp_ref in exp_refs:
                sample_ref = exp_to_sample.get(exp_ref)
                if sample_ref:
                    sample_to_files.setdefault(sample_ref, []).append(f)
            
            ana_refs = f.get("analysis", [])
            if isinstance(ana_refs, str):
                ana_refs = [a.strip() for a in ana_refs.split(",") if a.strip()]
            elif not isinstance(ana_refs, list):
                ana_refs = [ana_refs]
            
            for ana_ref in ana_refs:
                sample_to_files.setdefault(ana_ref, []).append(f)
    
    samples_data = []
    for sample in data.get("samples", []):
        sample_alias = sample.get("alias", "")
        sample_accession = sample.get("accession", "")
        
        ind_ref = sample.get("individual", "")
        individual_info = individuals.get(ind_ref, {})
        
        sample_id_clean = sanitize(sample_alias)
        individual_id_clean = sanitize(individual_info.get("alias", ind_ref))
        
        sex_raw = individual_info.get("sex")
        sex = sanitize(sex_raw) if sex_raw else "NA"
        
        disease_status_raw = str(sample.get("disease_or_healthy", "")).lower()
        status = 1 if "tumor" in disease_status_raw or "disease" in disease_status_raw else 0
        
        phenotype_raw = individual_info.get("phenotypic_features_terms", [])
        if isinstance(phenotype_raw, list):
            phenotype_str = ";".join([str(p) for p in phenotype_raw])
        else:
            phenotype_str = str(phenotype_raw)
        
        analysis_info = analyses.get(sample_alias, {})
        if not analysis_info and sample_accession:
            analysis_info = analyses.get(sample_accession, {})
            
        exp_method = sample_to_exp_method.get(sample_alias, "")
        if not exp_method and sample_accession:
            exp_method = sample_to_exp_method.get(sample_accession, "")
        
        files = []
        if sample_accession in sample_to_files:
            files.extend(sample_to_files[sample_accession])
        if sample_alias in sample_to_files:
            files.extend(sample_to_files[sample_alias])
        
        ana_files = analysis_info.get("research_data_files", [])
        if isinstance(ana_files, str):
            ana_files = [f.strip() for f in ana_files.split(",") if f.strip()]
        for f_ref in ana_files:
            if f_ref in data_files:
                files.append(data_files[f_ref])
            else:
                files.append({"name": f_ref})
        
        fastq_1, fastq_2 = [], []
        bams, bais = [], []
        crams, crais = [], []
        vcfs, others = [], []
        
        seen_files = set()
        for f in files:
            name = f.get("name", f.get("alias", ""))
            if not name or name in seen_files:
                continue
            seen_files.add(name)
            
            format_val = str(f.get("format", "")).upper()
            lower_name = name.lower()
            
            if input_directory:
                full_path = os.path.join(input_directory, name)
            else:
                full_path = name
            
            if format_val == "FASTQ" or lower_name.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fastq')):
                if '_r1' in lower_name or '_1.f' in lower_name:
                    fastq_1.append(full_path)
                elif '_r2' in lower_name or '_2.f' in lower_name:
                    fastq_2.append(full_path)
                else:
                    fastq_1.append(full_path)
            elif format_val == "BAM" or lower_name.endswith('.bam'):
                bams.append(full_path)
            elif format_val == "BAI" or lower_name.endswith('.bai'):
                bais.append(full_path)
            elif format_val == "CRAM" or lower_name.endswith('.cram'):
                crams.append(full_path)
            elif format_val == "CRAI" or lower_name.endswith('.crai'):
                crais.append(full_path)
            elif format_val == "VCF" or lower_name.endswith(('.vcf', '.vcf.gz', '.bcf')):
                vcfs.append(full_path)
            else:
                others.append(full_path)
        
        fastq_1.sort()
        fastq_2.sort()
        bams.sort()
        bais.sort()
        crams.sort()
        crais.sort()
        vcfs.sort()
        
        file_rows = []
        
        for f1, f2 in zip_longest(fastq_1, fastq_2):
            file_rows.append({
                "fastq_1": f1,
                "fastq_2": f2 if f2 else "",
                "single_end": "false" if f2 else "true"
            })
            
        for b, bi in zip_longest(bams, bais):
            file_rows.append({
                "bam": b,
                "bai": bi if bi else "",
                "single_end": "false"
            })
            
        for c, ci in zip_longest(crams, crais):
            file_rows.append({
                "cram": c,
                "crai": ci if ci else "",
                "single_end": "false"
            })
            
        for v in vcfs:
            file_rows.append({
                "vcf": v,
                "single_end": "false"
            })
            
        if others and not file_rows:
            file_rows.append({
                "data_files": ";".join(others),
                "single_end": "false"
            })
        
        if not file_rows:
            file_rows.append({
                "single_end": "false"
            })
        
        base_sample_data = {
            "sample": sample_id_clean,
            "individual_id": individual_id_clean,
            "sex": sex,
            "status": status,
            "phenotype": phenotype_str,
            "sample_type": sample.get("type", ""),
            "disease_status": sample.get("disease_or_healthy", ""),
            "case_control_status": sample.get("case_control_status", ""),
            "tissue": sample.get("biospecimen_tissue_term", ""),
            "experiment_method": exp_method,
            "analysis_method": analysis_info.get("analysis_method", "")
        }
        
        lane_counter = 1
        for r in file_rows:
            merged = {**base_sample_data, **r}
            merged["lane"] = "L" + str(lane_counter)
            lane_counter += 1
            samples_data.append(merged)

    all_columns = [
        "sample", "lane", "individual_id", "sex", "status", "phenotype", "sample_type",
        "disease_status", "case_control_status", "tissue", "experiment_method", "analysis_method",
        "fastq_1", "fastq_2", "single_end", "bam", "bai", "cram", "crai", "vcf", "data_files"
    ]
    
    samples_df = pd.DataFrame(samples_data)
    
    for col in all_columns:
        if col not in samples_df.columns:
            samples_df[col] = ""
            
    samples_df = samples_df[all_columns]
    
    samples_df.fillna("", inplace=True)
    samples_df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json", type=Path)
    parser.add_argument("--input_directory", type=str, default="")
    parser.add_argument("--output", type=str, default="samplesheet.csv")
    
    args = parser.parse_args()
    main(args.input_json, args.input_directory, args.output)