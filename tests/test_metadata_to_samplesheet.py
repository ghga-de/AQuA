import json
import tempfile
import unittest
from pathlib import Path

from bin.metadata_to_samplesheet import parse_metadata


class MetadataToSamplesheetTests(unittest.TestCase):
    def parse_inline_metadata(self, metadata):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.json"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
            return parse_metadata(str(metadata_path), "")

    def test_preserves_analysis_provenance_for_derived_files(self):
        rows = self.parse_inline_metadata(
            {
                "individuals": [{"alias": "ind_1", "sex": "FEMALE"}],
                "samples": [
                    {
                        "alias": "sample_1",
                        "individual": "ind_1",
                        "disease_or_healthy": "Healthy",
                    }
                ],
                "experiment_methods": [
                    {
                        "alias": "em_wgs",
                        "instrument_model": "ILLUMINA",
                        "library_type": "WGS",
                        "sequencing_layout": "PE",
                    }
                ],
                "analysis_methods": [
                    {
                        "alias": "am_align",
                        "type": "WGS",
                        "software": "bwa-mem2",
                        "reference_genome": "GRCh38",
                        "target_bed": "capture.bed",
                    }
                ],
                "experiments": [
                    {
                        "alias": "exp_1",
                        "sample": "sample_1",
                        "experiment_method": "em_wgs",
                    }
                ],
                "research_data_files": [
                    {
                        "alias": "r1",
                        "name": "sample_1_R1.fastq.gz",
                        "format": "FASTQ",
                        "experiments": ["exp_1"],
                    }
                ],
                "analyses": [
                    {
                        "alias": "analysis_1",
                        "analysis_method": "am_align",
                        "research_data_files": ["r1"],
                        "analysis_tool": "sentieon",
                        "analysis_genome": "GRCh37",
                        "alignment_status": "unaligned",
                    }
                ],
                "process_data_files": [
                    {
                        "alias": "bam",
                        "name": "sample_1.bam",
                        "format": "BAM",
                        "analysis": "analysis_1",
                    }
                ],
            }
        )

        bam_row = next(row for row in rows if row.get("bam"))
        self.assertEqual(bam_row["analysis_method"], "wgs")
        self.assertEqual(bam_row["analysis_tool"], "sentieon")
        self.assertEqual(bam_row["analysis_genome"], "GRCh37")
        self.assertEqual(bam_row["target_bed"], "capture.bed")
        self.assertEqual(bam_row["alignment_status"], "unaligned")

    def test_pairs_fastqs_with_underscore_one_two_suffixes(self):
        rows = self.parse_inline_metadata(
            {
                "individuals": [{"alias": "ind_1"}],
                "samples": [{"alias": "sample_1", "individual": "ind_1"}],
                "experiment_methods": [
                    {
                        "alias": "em_wgs",
                        "instrument_model": "ILLUMINA",
                        "library_type": "WGS",
                        "sequencing_layout": "PE",
                    }
                ],
                "experiments": [
                    {
                        "alias": "exp_1",
                        "sample": "sample_1",
                        "experiment_method": "em_wgs",
                    }
                ],
                "research_data_files": [
                    {
                        "alias": "r2",
                        "name": "sample_1_2.fastq.gz",
                        "format": "FASTQ",
                        "experiments": ["exp_1"],
                    },
                    {
                        "alias": "r1",
                        "name": "sample_1_1.fastq.gz",
                        "format": "FASTQ",
                        "experiments": ["exp_1"],
                    },
                ],
            }
        )

        fastq_row = next(row for row in rows if row.get("fastq_1"))
        self.assertEqual(fastq_row["fastq_1"], "sample_1_1.fastq.gz")
        self.assertEqual(fastq_row["fastq_2"], "sample_1_2.fastq.gz")
        self.assertEqual(fastq_row["single_end"], "false")


if __name__ == "__main__":
    unittest.main()
