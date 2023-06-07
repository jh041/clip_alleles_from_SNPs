# clip_alleles_from_SNPs
Noisy alleles can get your SNP filtered from the data. Clip them first. Save the SNP!

This script is for an R function that reads in a .vcf file of genotypes, calculates allele frequencies, identifies low-frequency alleles, and removes the allele from the site without deleting the allele. Genotypes of the low-frequency alleles are replaced with missing data, and a new .vcf is outputed.

Imagine you have a trinary SNP with allele frequencies like this...

G: 0.8704, C: 0.1284, A: 0.0012

This might be an important SNP, but the site won't survive genotyping filtering if you are excluding sites that are non-binary SNPs, or have minor allele frequencies less than 0.01.

However, after allele clipping that same site will pass most filters...

G: 0.871291,	C: 0.128709

The function is written mostly in base R, and only requires the vcfR package to parse the .vcf file.

Usage is straightforward:

```
clip_alleles("input_vcf_file_name.vcf", "output_vcf_file_name.vcf.gz", maximum frequency of noisy alleles, e.g. 0.01)
```

A few warnings:

This function presently only works on unphased .vcf genotypes. Phased genotypes will cause the function to throw an error.

This function only clips ALT alleles with low frequencies from sites. If the frequency of the REF allele falls below the max allele frequency, it is replaced in the genotypes by missing data but not clipped. Doing that would conflict with the genomic reference sequence and possibly cause issues in downstream applications like microhaplotyping.

While we're talking about this, if your genomic reference comes from the same population as your samples, and the reference allele
isn't present in the population genotypes, then you probably have a sequencing error in your reference.

Subsequent maf filtering will remove garbage sites after processing, including those that are left monomorphic. Therefore, it is recommended that you set the max allele frequency variable no higher than what you are willing to accept as a minor allele frequency for your data.

So far, this script has been successfully tested on vcf outputs from three SNP callers: FreeBayes, GATK Haplotype caller, and BCFtools, but should work on most vcf files, as long as the genotypes are unphased.

Also, missing-data filters used afterwards could exclude these loci of interest if the max allele frequency is set too high.
