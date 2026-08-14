#!/usr/bin/env python3
import csv
import sys
import os

samplesheet_path = sys.argv[1]
unified_all_path = sys.argv[2]
multiqc_data_dir = sys.argv[3]

# 1. Parse Tool Versions Flexibly
tool_versions = {}
versions_file = os.path.join(multiqc_data_dir, 'multiqc_software_versions.txt')
if os.path.exists(versions_file):
    with open(versions_file, mode='r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = []
        for idx, row in enumerate(reader):
            if idx == 0:
                headers = [h.strip().lower() for h in row]
            else:
                for col_idx in range(1, len(headers)):
                    if col_idx < len(row) and row[col_idx].strip():
                        tool_versions[headers[col_idx]] = row[col_idx].strip()
                        tool_versions[row[0].strip().lower()] = row[col_idx].strip()

# 2. Parse Samplesheet Dynamically and Merge Duplicate Sample Rows
merged_samples = {}
with open(samplesheet_path, mode='r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        row_data = {k.strip().lower(): v.strip() for k, v in row.items() if k}
        s = row_data.get('sample', '')
        if s:
            if s not in merged_samples:
                merged_samples[s] = row_data
            else:
                # If sample exists, update empty fields with new files (e.g. adding fastq to a bam row)
                for k, v in row_data.items():
                    if v and not merged_samples[s].get(k):
                        merged_samples[s][k] = v

samples_data = list(merged_samples.values())

def get_file_for_sample(row_data, mqc_samp, tool_name):
    t_lower = tool_name.lower()
    
    # Determine the priority file types based on the tool
    if 'fastp' in t_lower or 'nanostat' in t_lower or 'fastqc' in t_lower:
        priorities = ['fastq_1', 'fastq_2']
    elif 'mosdepth' in t_lower or 'samtools' in t_lower or 'cramino' in t_lower:
        priorities = ['bam', 'cram']
    elif 'bcftools' in t_lower:
        priorities = ['vcf']
    else:
        priorities = ['fastq_1', 'bam', 'cram', 'vcf']

    # Step A: Check priority files to see if the MultiQC name matches the filename
    for key in priorities:
        val = row_data.get(key, '')
        if val and (mqc_samp in val or val in mqc_samp):
            return val
            
    # Step B: If no string match, return the first available priority file for this tool
    for key in priorities:
        val = row_data.get(key, '')
        if val:
            return val
            
    # Step C: Absolute fallback if tool priorities fail
    for key in ['fastq_1', 'bam', 'cram', 'vcf']:
        if row_data.get(key):
            return row_data[key]
            
    return ""

# 3. Process the MultiQC Data
with open(unified_all_path, mode='r', encoding='utf-8') as fin, \
     open('qc_metrics_unified.tsv', mode='w', encoding='utf-8', newline='') as fout:
    
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(['sample', 'file', 'tool', 'tool_version', 'qc_metrics', 'definition', 'value'])

    reader = csv.DictReader(fin, delimiter='\t')
    
    for row in reader:
        concept = row.get('concept', '').strip()
        if not concept:
            continue
            
        mqc_samp = row.get('sample', '').strip()
        tool = row.get('tool', '').strip()
        metric = row.get('metric', '').strip()
        val_raw = row.get('val_raw', '').strip()
        description = row.get('description', '').strip()

        input_file = ""
        real_sample_name = mqc_samp
        
        # Match Sample ID to File path
        for s in samples_data:
            if s['sample'] == mqc_samp:
                real_sample_name = s['sample']
                input_file = get_file_for_sample(s, mqc_samp, tool)
                break
                
        if not input_file:
            sorted_samples = sorted(samples_data, key=lambda x: len(x['sample']), reverse=True)
            for s in sorted_samples:
                if s['sample'] in mqc_samp or mqc_samp in s['sample']:
                    real_sample_name = s['sample']
                    input_file = get_file_for_sample(s, mqc_samp, tool)
                    break
                    
                for key in ['fastq_1', 'fastq_2', 'bam', 'cram', 'vcf']:
                    val = s.get(key, '')
                    if val and (mqc_samp in val or val in mqc_samp):
                        real_sample_name = s['sample']
                        input_file = val
                        break
                if input_file:
                    break

        # Match Tool Version
        t_lower = tool.lower()
        base_tool = t_lower.split('_')[0]
        version = tool_versions.get(t_lower, "")
        
        if not version:
            version = tool_versions.get(base_tool, "")
            
        if not version:
            for v_tool, v_ver in tool_versions.items():
                if t_lower in v_tool or v_tool in t_lower or base_tool in v_tool:
                    version = v_ver
                    break

        writer.writerow([real_sample_name, input_file, tool, version, metric, description, val_raw])