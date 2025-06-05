VCF Filtering and Annotation Pipeline
A comprehensive Python-based pipeline for filtering somatic variants from tumor-only VCF files and performing clinical annotation with gene and disease association analysis.

Overview
This pipeline consists of two main components:

VCF Filtering (vcf_filter2.py): Applies multi-layered quality filters optimized for tumor-only somatic variant detection
Clinical Annotation (clinical_test_annotation.py): Annotates filtered variants with gene information, clinical significance, and disease associations
The filtering strategy is specifically designed for tumor-only samples, addressing the primary challenge of distinguishing somatic variants from germline variants without matched normal tissue.

# Requirements
1. System Requirements
Python 3.7 or higher
2. Python Dependencies
pip install pysam
3. Input File Requirements
; VCF File Format
Your input VCF file must contain the following INFO fields for filtering:
TLOD: Tumor Log Odds score (Float)
DP: Total read depth (Integer)
POPAF: Population allele frequency on -log10 scale (Float)
GERMQ: Phred-scaled germline confidence (Integer)
AF: Allele frequency in tumor (Float)
CONTQ: Contamination quality score (Integer)
SEQQ: Sequencing quality score (Integer)
MPOS: Median position in reads (Integer)
; JSON Criteria File
The multiallelic_criteria.json file defines filtering thresholds and is included in the pipeline.

# Usage
1. Filter Raw VCF File
python scripts/vcf_filter2.py -i input.vcf -o output_filtered.vcf -c multiallelic_criteria.json
; With Verbose Output
python scripts/vcf_filter2.py -i input.vcf -o output_filtered.vcf -c multiallelic_criteria.json --verbose

Command Line Arguments
-i, --input: Input VCF file path (required)
-o, --output: Output filtered VCF file path (required)
-c, --criteria: JSON criteria file path (required)
-v, --verbose: Enable verbose logging (optional)

Example Commands
; Filter Mutect2 Raw Output:
python scripts/vcf_filter2.py -i SAMPLE_mutect2_raw.vcf -o mutect2_filtered.vcf -c multiallelic_criteria.json --verbose

; Filter Clinical Test Data:
python scripts/vcf_filter2.py -i clinical_test.vcf -o clinical_test_filtered.vcf -c multiallelic_criteria.json --verbose

2. Annotate Filtered Variants
python clinical_test_annotation.py
Note: The annotation script is currently configured to process clinical_test_filtered.vcf by default. To process different files, modify the vcf_file variable in the main() function.

; Annotation Features: 
Gene identification using genomic coordinates
Clinical significance assessment
Disease association lookup
Allele frequency annotation
Excel and CSV output formats
Filtering Criteria Details
High Confidence Tumor-Only Strategy


; Filter	Threshold	Purpose:
TLOD >=20.0	Statistical confidence that variant exists in tumor
DP	>=50	Minimum read depth for reliable variant calling
POPAF	>=3.0	Population frequency filter (0.1% threshold on -log10 scale)
GERMQ	>=30	99.9% confidence variant is somatic, not germline
AF	>=0.05	Minimum 5% allele frequency to exclude artifacts
CONTQ	>=20	99% confidence variant is not contamination
SEQQ	>=20	99% confidence variant is not sequencing error
MPOS	>=10	Excludes variants within 10bp of read ends
（PASS variants meet all criteria）

; Expected Filtering Outcomes
Germline variant removal: 95-98%
Artifact reduction: 85-90%
True somatic retention: 85-92%
Manual review candidates: 2-5% of filtered variants
Overall retention: 8-15% of raw variants

; Filter Hierarchy
Tier 1 (Critical): POPAF, GERMQ - Address fundamental tumor-only challenges
Tier 2 (Quality): TLOD, DP, SEQQ - Ensure statistical confidence
Tier 3 (Artifact): CONTQ, MPOS, AF - Remove technical artifacts

# Testing
; Test Data Preparation
1. Create Test VCF File
The pipeline includes clinical_test.vcf with 10 test variants representing different scenarios:

# Verify test file format
head -20 clinical_test.vcf

2. Validate JSON Criteria
# Check JSON syntax
python -m json.tool multiallelic_criteria.json

; Running Tests
1. Test VCF Filtering
# Run filtering with verbose output
python scripts/vcf_filter2.py -i clinical_test.vcf -o clinical_test_filtered.vcf -c multiallelic_criteria.json --verbose

Expected Output:
2024-XX-XX XX:XX:XX - INFO - Processed 10 records, X passed.

2. Test Clinical Annotation
# Run annotation (requires internet connection)
python clinical_test_annotation.py
Expected Output:

Analyzing VCF file: clinical_test_filtered.vcf
============================================================
Analyzing variant: 12:25398284 C>T
Analyzing variant: 17:7674220 C>A
...
Found X PASS variants
Processing results...
Excel file saved: clinical_test_annotation_results.xlsx
CSV backup saved: clinical_test_annotation_results.csv

; Output Files Generated:
clinical_test_annotation_results.xlsx - Formatted Excel report
clinical_test_annotation_results.csv - CSV backup

; Expected Annotation Fields:
Chromosome, Position, RS_ID
Reference/Alternate alleles
Gene name, ID, and description
Allele frequency
Clinical significance
Associated diseases

； The annotation script uses the following external APIs:
Ensembl REST API: Gene and variant information
NCBI E-utilities: ClinVar disease associations


# Troubleshooting Tests
Common Issues and Solutions
1. Missing Dependencies
# Error: No module named 'pysam'
pip install pysam
# Error: No module named 'openpyxl'
pip install openpyxl  # Optional for Excel output

2. Network Issues in Annotation
# If APIs are unreachable, annotation will show:
Error fetching gene for X:Y - Connection timeout
Solution: Check internet connection and retry

3. VCF Format Issues
# Error: Invalid condition for field 'X'
Solution: Verify VCF contains all required INFO fields

4. Empty Results
# If no variants pass filtering:
Processed X records, 0 passed.
Solution: Review filtering criteria or input data quality


# Limitations
Designed specifically for tumor-only samples
Requires specific INFO field annotations
Internet connection required for clinical annotation
API rate limiting may slow annotation of large datasets
Clinical significance interpretations should be validated independently
