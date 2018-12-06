# PhasedPileup
Make a display with reads from a BAM file phased to the correct parent genome
## Example

![Example output](doc/General.png)

## Requirements and Usage

This package requires:
* pysam  (tested with 0.12.0.1)
* svgwrite (tested with 1.1.11)

The most common usage is to view the pileup of a particular gene:

    # Phased reads around a specified gene (YFG01 = your favorite gene 1)
    python PhasedPileup.py --gtf-file ANNOTATION.gtf --gene YFG01 snp_file.txt samfile.bam
which will then output `YFG01_samfile_phased.svg` (though there are options to
control this).

But you can also provide any arbitrary genomic region:

    # Phased reads in a genomic region
    python PhasedPileup.py --coords chr1:3600-3800 snp_file.bed samfile.bam
    python PhasedPileup.py --coords chr9:3600..3800 snp_file.bed samfile.bam

The SNP file is assumed to be in the format:

    CHROM	POS0 POS1	REF|ALT

where POS0 is the 0-based coordinate and POS1 is the 1-based coordinate of the
SNP. It is furthermore assumed that there are no indels in the file... I
haven't considered what might happen if they were in there.


## Colors and visual style

Red reads are the reference allele, and blue reads are the alternate allele.
Gray reads are reads that cannot be unambiguously phased, either due to not
overlapping a SNP, or having SNPs with conflicting or unknown phase.

![Ref/Alt Example](doc/RefAlt.png)

Additionally, read1 is marked with a black flag at the start of the read.

![Read 1 is marked with a flag](doc/Read1Flag.png)


Reads that span an intron (or, more generally, with a gap relative to the
reference sequence) are shown with a horizontal line in the phase color that
spans the intron. If you have paired-end reads, the unsequenced portion of the
insert is shown with a black line connecting the tips of the two reads.  If the
forward and reverse read overlap, then the black line will instead indicate the
overlapping portion of the read.

At the moment, these are all hard-coded in, but they can either be edited in
the code or by editing the CSS style information of the produced SVG file.

