"""
sam_to_bed.py

Searches through 1) a file of SNP reads and 2) a gzipped SAM file of RNAseq
reads to output a gzipped BED file with reads colored according to the SNPS
within them. SNP locations within the reads will also be output in Uppercase.

Reads will be colored...
Red if SNP in BL6
Blue if SNP in Cast
Black if SNP that exists in neither BL6 nor Cast is found, or if no SNPs at all

"""
import sys
import gzip

samfile, snpfile = sys.argv[1:]

snpdic = {}
# Initialize Snp counting var and other vars needed for printing.

cast_track = 0
bl_track = 0
neither_track = 0
both_track = 0
bed_seq = ""
bed_color = "0,0,0"
bed_strand = ""

# Code constants into variables?

print("track name=F1_RNA-seq visibility=squish itemRgb=On db=mm9") # Why is this constant in code bad?

with open(snpfile, mode="rt") as snps:
    for snpline in snps:
        snp_chrom, snp_pos, genes, bl, cast, prop, pile, basequal \
            = snpline.rstrip("\r\n").split("\t")
        bl = bl.upper()
        cast = cast.upper()
        snp_pos = int(snp_pos)
        snp_chrom = "chr" + snp_chrom

        # Creating dictionary of all required information mapped to the posi-
        # tion of the key.
        snpdic[snp_pos] = (snp_chrom, bl, cast)

with gzip.open(samfile, mode="rt", encoding="utf-8") as sams:
    for line in sams:
        a, strand, samchrom, start, c, d, e, f, g, seq, h, i, j, k \
            = line.rstrip("\r\n").split("\t")

        # Putting into correct forms. The seq is all lowercase as non-SNP bases
        # will be in lowercase printed later on.
        seq = seq.lower()
        start = int(start)

        for i, base in enumerate(seq):

            # Allows us to keep track of the location of SNP
            basepos = start + i

            if basepos in snpdic.keys():

                # Check that the chromosomes match between sam and snp files.
                # If not, string will not be included in the output file.
                if snpdic[basepos][0] != samchrom:
                    break

                bl_test = snpdic[basepos][1]
                cast_test = snpdic[basepos][2]

                if base == bl_test.lower():
                    bed_seq += bl_test
                    bl_track += 1
                elif base == cast_test.lower():
                    bed_seq += cast_test
                    cast_track += 1
                else:
                    # Makes sure that the non-ref(Bl) or cast SNPs are also
                    # uppercase.
                    bed_seq += base.upper()
                    neither_track += 1
            else:
                bed_seq += base

        # Color matching. Each statement in the if statements will point to a
        # certain condition
        if bl_track == 0 and cast_track == 0 and neither_track != 0:
            bed_color = "0,0,0"
        elif bl_track != 0 and cast_track == 0 and neither_track == 0:
            bed_color = "255,0,0"
        elif bl_track == 0 and cast_track != 0 and neither_track == 0:
            bed_color = "0,0,255"

        # Assign strand fwd/reverse depending on given variables
        if strand == "0":
            bed_strand = "+"
        elif strand == "16":
            bed_strand = "-"

        # Correcting for the fact that bed files are 0-indexed, not 1-indexed
        bed_start = start - 1
        bed_stop = bed_start + len(seq)

        print(samchrom, bed_start, bed_stop, bed_seq, "0", bed_strand,
              "0", "0", bed_color, sep="\t")
