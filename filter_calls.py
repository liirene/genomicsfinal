"""
filter_calls.py

Applies filters to calls to get rid of 1) sequencing biases 2) calling biases
3) any other biases that I manage to get out

input: VCF file with calls?
output: 1) VCF file with calls filtered 2) Stats file with reasons for each of
the ones that were taken out to be taken out lel

"""

# Something here to decide whether running as a module or as a main system file
import sys
import re

# Currently call looks like: python3 filtercalls.py input.vcf output.vcf

# define variables
inputvcf = sys.argv[1]

# empty dictionary for calls and their associated attributes
calls = {}

# filter out all calls that don't have a call quality associated
with open(inputvcf, mode = "rt") as in_vcf:
    for line in in_vcf:
        if not line.startswith("#"):
            chrom, pos, ntid, ref, alt, qual, filt, info, line_format, \
            var_match = line.rstrip(
                "\t\n").split("\t")

            # Only focus on SNP values
            if not qual == ".":
                info_dict = {}
                clinvar_query = []

                # pull out query values for omim

                # Pull out query values for clinvar
                num = re.findall("(?<=chr)\w", chrom)[0] + "[chr]"
                gene = re.findall("\(([^\)]+)\)", var_match)
                hgnv = re.findall("c\.[^\,\:\(\s]+", var_match)
                refseq = re.findall("N\w?_[^\(]+", var_match)

                query_list = [gene, hgnv, refseq]

                # test if any are none

                clinvar_query.append(num)
                for item in query_list:
                    if not len(item) == 0:
                        clinvar_query.append(item[0])

                # Archive the information for use later... make into lists
                info = info.split(";")
                line_format = line_format.split(":")
                var_match = var_match.split(":")

                clinvar_query = " AND ".join(clinvar_query)

                # Create dictionary with important values for our analyses
                for i, field in enumerate(line_format):
                    info_dict[field] = var_match[i]

                print(clinvar_query)

                #calls.update(
                  #  {(chrom, pos): [ref, alt, qual, info_dict, clinvar_query]})
