clip_alleles <- function(path_2_vcf_file, path_2_new_vcf_file.gz, max_allele_frequency) {

require(vcfR)

### Read in vcf data
vcf <- read.vcfR(path_2_vcf_file, verbose=F)

### extract genotypes
genotypes <- extract.gt(vcf)

### You now have an object with genotypes that looks like this:
# "0/0" "1/1"  "2/2" "0/2"
# "0/1" "1/1" "1/2" NA
### This is easier to work with than the
### way genotypes appear in the raw vcfR object
# > vcf@gt[1:3,1:2]
# "GT:AD:DP:GQ:PL" "1/1:0,1,0:1:3:41,3,0,41,3,41"
# "GT:AD:DP:GQ:PL" "./.:.:.:.:."
### So we'll use the simple matrix to inform changes
### to the more complex one that will eventually be 
### written to file as the new .vcf

### Note that the first column in vcf@gt is not a genotype
### So the coordinates between "genotypes" and vcf@gt are off
### We'll add an extra column at the front of "genotypes"
### To correct for this.

genotypes <- cbind(vcf@gt[,1],genotypes)

### Extract what a missing entry looks like 
### for this particular vcf file 

missing <- vcf@gt[is.na(genotypes)][1]

### Extract a list of the Reference and alternate allelels

nucleotides <- getFIX(vcf)[,4:5]

### Loop through sites, calculate an allele frequency table for each locus
#       A 0 0.3709524
#       T 1 0.2495238
#       G 2 0.3795238
### 1) Replace the genotype in vcf@gt as missing if one allele is less than "max_allele_frequency"
### 2) Renumber the alleles in genotypes, if the clippled allele is not the last alternate allele
### 3) Remove clipped alleles from the ALT columns of the vcfR object

for(site in 1:length(row.names(genotypes))) {

   locus <- unname(genotypes[site,-1])
   ref.alt <- c(unname(nucleotides[site,1]), unname(unlist(strsplit(nucleotides[site,2],","))))
   variants <- unname(unlist(sapply(locus, function(x) strsplit(x, "/"))))
   new_record <- data.frame(ref.alt, 0:(length(ref.alt) - 1), rep(0, length(ref.alt)))
   
   for(record in 1:length(new_record[,2])){
   
      freq <- table(variants[which(variants == new_record[record, 2])])/sum(table(variants))
      if(length(freq) != 0) { freq -> new_record[record,3] }
   }
   
   ### Step 1
   low_freq_alleles <- new_record[which(new_record[,3] < max_allele_frequency), 2]
   
   for(allele in low_freq_alleles) {
   
      vcf@gt[site, grep(allele, genotypes[site,])] <- missing
   
   }
   
   #Step 2
   high_freq_nucleotides <- new_record[which(new_record[,3] > max_allele_frequency), 1]
   
   if(length(low_freq_alleles) != 0 & high_freq_nucleotides[1] == ref.alt[1]) {
   
      old_alt_allele_numbers <- new_record[which(new_record[,2] > max(low_freq_alleles)),2]
   
      if(!(length(old_alt_allele_numbers) == 0 & low_freq_alleles[1] == 0)) {
   
         for(number in old_alt_allele_numbers) {
      
            vcf@gt[site, grep(paste0(number,"/"),vcf@gt[site,])] <- gsub(paste0(number,"/"), paste0(number - 1,"/"), vcf@gt[site, grep(paste0(number,"/"),vcf@gt[site,])])
            vcf@gt[site, grep(paste0("/",number),vcf@gt[site,])] <- gsub(paste0("/", number), paste0("/", number - 1), vcf@gt[site, grep(paste0("/", number),vcf@gt[site,])])
      
         }
   
      }
   }
   
   #Step 3   
   if(high_freq_nucleotides[1] != ref.alt[1]) {
   
      next   # The reference allele is absent in the genotypes. Don't replace with an alt allele.
   
   }
   
   else if(length(high_freq_nucleotides) > 1) {
   
      vcf@fix[site, 4] <- high_freq_nucleotides[1]
      vcf@fix[site, 5] <- paste(high_freq_nucleotides[2:length(high_freq_nucleotides)], collapse=",")
   
   }
}

### write a new vcf to file

write.vcf(vcf, file = path_2_new_vcf_file.gz)

}
