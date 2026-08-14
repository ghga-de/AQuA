#!/usr/bin/env python3
import csv
import sys

samplesheet_path = sys.argv[1]
unified_all_path = sys.argv[2]

samples_data = {}
with open(samplesheet_path, mode='r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        s_id = row.get('sample', '').strip()
        if s_id:
            samples_data[s_id] = {
                'fastq_1': row.get('fastq_1', '').strip(),
                'bam': row.get('bam', '').strip(),
                'cram': row.get('cram', '').strip(),
                'vcf': row.get('vcf', '').strip()
            }

with open(unified_all_path, mode='r', encoding='utf-8') as fin, open('qc_metrics_unified.tsv', mode='w', encoding='utf-8', newline='') as fout:
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(['sample', 'file', 'tool', 'qc_metrics', 'definition', 'value'])
    
    reader = csv.reader(fin, delimiter='\t')
    next(reader)
    
    for row in reader:
        if len(row) > 1 and row[1]:
            mqc_samp = row[0].strip()
            tool = row[2].strip()
            metric = row[5].strip()
            val = row[6].strip()
            defn = row[10].strip() if len(row) > 10 else ""

            input_file = ""
            
            # 1. Try exact match first
            if mqc_samp in samples_data:
                s_info = samples_data[mqc_samp]
                for key in ['bam', 'cram', 'fastq_1', 'vcf']:
                    if s_info.get(key):
                        input_file = s_info[key]
                        break
            
            # 2. Fallback: match if sample name appears anywhere in the file path
            if not input_file:
                for s_id, s_info in samples_data.items():
                    for key, f_path in s_info.items():
                        if f_path and (mqc_samp in f_path or s_id in mqc_samp or mqc_samp in s_id):
                            input_file = f_path
                            break
                    if input_file:
                        break

            writer.writerow([mqc_samp, input_file, tool, metric, defn, val])