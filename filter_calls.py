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
