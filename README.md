# Description
The Barcode-Split algorithm is encapsulated within the "split_bam_by_tag" function of the Python package bebam.

# Dependencies for the bebam Package:

1. Python version 3.6 or higher
2. Samtools (must support the "samtools view –tag-file" and "samtools sort -t" parameters)
3. Bgzip
4. Mawk
5. Tabix


# Installation of the bebam Package

1. Download the package using "git clone https://github.com/linguowang007/Barcode-Split.git".
2. Navigate to the bebam folder.
3. Important!!! Modify lines 9-12 in the split_bam_file.py script within the bebam folder to replace the values of the samtools, bgzip, mawk, and tabix variables with the absolute paths of these tools in your operating system.
4. Return to the parent directory and execute "python setup.py install" in the terminal to complete the installation of the bebam package.


# How to Use Barcode-Split

The Barcode-Split algorithm, encapsulated in the split_bam_by_tag function within the bebam package, can be utilized in two ways:

## As a Standalone Command-Line Tool. 
After installing bebam, the split_bam_by_tag tool is automatically added to your system's environment variables. Use "which split_bam_by_tag" to locate it.

### Usage:
split_bam_by_tag [-h] --bam BAM --tag_list TAG_LIST --out_dir OUT_DIR [--nt NT] [--tag TAG] [--tag_type TAG_TYPE]

### Example command in the terminal:
split_bam_by_tag --bam bamfile.bam --tag_list tags.txt --out_dir output-dir-name --nt 8 --tag CB --tag_type Z

This command will split bamfile.bam by the values of CB and output the sub-BAM files to output-dir-name.

## As a Python Package Called from Other Python Scripts. 
Add the following code to your Python script (or an interactive Python session):
### test.py：
from bebam import split_bam_by_tag

help(split_bam_by_tag) # help infromation of split_bam_by_tag function

split_result = split_bam_by_tag(bam="path/to/bamfile.bam",
                 tag_list="path/to/tags.txt",
                 out_dir="out-dir-name",
                 nt=8,
                 tag="CB",
                 tag_type="Z")
                 
#"split_result" is a return dictionary where keys are CBs and values are paths to the split sub-BAM files.


# Run test on Demuxlet paper data

wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam

samtools index -@ 8 A.merged.bam

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz

gunzip GSM2560245_barcodes.tsv.gz

split_bam_by_tag --bam A.merged.bam --tag_list GSM2560245_barcodes.tsv --out_dir split-test
