"""
create_visualfiles.py

Creates a bigBED file capable of being loaded into any conventional genome
viewer.

Colors are given based on...
     Type of mutation

Input: file of SNP / mut reads and file of SAM/BAM reads
Output: bigBED file with fields colored

"""
# Include information in header about why the calls are colored as they are

# Keep in mind the indexing of sam/bam files!! Modify sam_to_bed.py for this
# purpose

# create bed detail file
# Header text: 'track name = <name> type=bedDetail description="desc" usescore=1 db=<hg19> visibility=3'
