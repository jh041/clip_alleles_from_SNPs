# clip_alleles_from_SNPs
Noisy alleles can get your SNP filtered from the data. Clip them first. Save the SNP!

This script is for an R function that reads in a .vcf file, calculates allele frequencies from genotypes, identifies low-frequency alleles, and removes the allele from the SNP or indel without deleting the site. Genotypes of the low-frequency alleles are replaced with missing data, and a new .vcf is outputed.

Imagine you have a trinary SNP with allele frequencies like this...

G: 0.8704, C: 0.1284, A: 0.0012

This might be an important SNP, but the site won't survive genotype filtering if you are excluding sites that are non-binary variants, or have minor allele frequencies less than 0.01.

However, after allele clipping that same site will become a binary SNP that will pass common minor allele frequency criteria...

G: 0.871291,	C: 0.128709

The function is written mostly in base R, and only requires the vcfR package to parse the .vcf file.

Usage is straightforward:

```
clip_alleles("input_vcf_file_name.vcf", "output_vcf_file_name.vcf.gz", maximum frequency of noisy alleles, e.g. 0.01)
```

A few considerations:

This function presently only works on unphased .vcf genotypes. Phased genotypes will cause the function to throw an error.

This function only clips ALT alleles with low frequencies from sites. If the frequency of the REF allele falls below the max allele frequency, it is replaced in the genotypes by missing data but not clipped. Doing that would conflict with the genomic reference sequence and possibly cause issues in downstream applications and interpretations of results.

No sites are deleted by the function. Removing sites afterwards must be done with additional filtration steps, such as minor allele frequency filtering. Therefore it is recommended that the maximum allele frequency variable be set no higher than what you are willing to accept as a minor allele frequency.

So far, this script has been successfully tested on vcf outputs from three SNP callers: FreeBayes, GATK Haplotype caller, and BCFtools, but should work on most vcf files, as long as the genotypes are unphased.

Also, be warned that missing-data filters used afterwards could exclude these loci of interest if the max allele frequency is set too high.
