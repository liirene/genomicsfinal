"""
snp_ids.py

Queries databases to find information for SNPs

input: VCF? BAM? file with SNP information (in terms of location)

output: VCF file with SNPs annotated with gene, probable effect, etc.

** Make sure to put organism here and specify the build (through input?)

"""

import sys
import urllib2
