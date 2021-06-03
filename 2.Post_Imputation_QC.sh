### Post-imputation QC

# Download the data from the server 
# To unzip the data:
unzip -P PASSWORD '*.zip'

# Generating Softcall + hardcall binaries:
# Softcall: rsq > 0.3
# Hardcall: rsq > 0.8

# Optional Step: what sites pass our MAF and Rsq cutoffs? (Only to check)
Copy and save the following as filter_variants_maf_r2_updated.R

if (!"tidyverse" %in% rownames(installed.packages())) {
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
}
if (!"data.table" %in% rownames(installed.packages())) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(data.table)

arg <- commandArgs(trailingOnly = T)

# define function we will use
filter_info <- function (chr, maf, rsq) {
  # Read the info file of specific chromosome
  data <- paste0("chr", chr,".info") %>% fread()
  # Filter it to maf and rsq
  dat <- data[data$MAF >= maf & data$Rsq >= rsq]
  # Generate chromosome, bp, and range from the "SNP" column
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste0(dat$chr, ":", dat$bp, "-", dat$bp)
  # writing files
  dat[,c("SNP","ALT_Frq","Rsq")] %>% fwrite(
    paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".info"),
    row.names = F, quote = F, sep = "\t"
  )
  dat[,c("range")] %>% fwrite(
    paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".txt"),
    col.names = F, row.names = F, quote = F, sep = "\t"
  )
  # report number of SNPs that pass the filter
  paste0("Number of SNPs in chromosome ", chr, " after filter (MAF >= " , maf, ", Rsq >= ", rsq, "): ", nrow(dat)) %>% print()
}

# SET MAF and RSQ filters here:

lapply(1:22, filter_info, maf = arg[1], rsq = arg[2])

gunzip *.info.gz

# SYNTAX: 
Rscript filter_variants_maf_r2_updated.R [MAF cutoff] [Rsq cutoff]
# following will use MAF cutoff of 0.05 and Rsq cutoff of 0.3:
# Rscript filter_variants_maf_r2_updated.R 0.05 0.3


# Now let's generate a softcall. The script below converts the dose.vcf.gz files to Plink binaries while filtering by given Rsq cutoff value.

# make the output folder
mkdir plink_files_soft (only rsq filter, no man filter here)

# change "--qual-threshold" to desired rsq threshold. in this case, it is 0.3 because we are generating softcall

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink --vcf chr$chnum.dose.vcf.gz --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.3 --make-bed --out plink_files_soft/chr$chnum  --double-id
done

cd plink_files_soft
ls | grep ".bim" > merge_list.txt
# sed -i 's/.bim//g' merge_list.txt Get "sed: 1: "merge_list.txt": invalid command code m", using the code below instead
sed -i '' -e 's/.bim//g' merge_list.txt 
plink --merge-list merge_list.txt --make-bed --out MPBC.Rsq03

#NOTE! No phenotypes, sex or MAF filtering in this file! --> Add pheno, sex and set MAF threshold before running GWAS
#The IID and FID is different in the pre-imputed file and post-imputed file, this need to be changed in the pre-imputed fam file (create a new) in order to update the info.

# In R:
## Creating new file for updating the pheno in the PLINK files (.bed, .bim) after imputation (the ID changed for me after imputation, FID and IID merged)
setwd("~/Desktop/MPBC/plink_files_soft")
# fam file after imputation with no phenotypes (to check)
MPBC.Rsq03 <- read.table("~/Desktop/MPBC/plink_files_soft/MPBC.Rsq03.fam", quote="\"", comment.char="")
# fam file before imputation with correct phenotypes etc
FILTERED.MPBC <- read.table("~/Desktop/MPBC/FILTERED.MPBC.fam", quote="\"", comment.char="")
# Create a column with correct IID and make a new fam file
FILTERED.MPBC$IID_new <- paste0(FILTERED.MPBC$V1, "_", FILTERED.MPBC$V2)
UpdatedIDs <- FILTERED.MPBC[,c("IID_new","IID_new", "V3", "V4", "V5", "V6")]
write.table(UpdatedIDs, file=paste("UpdatedIDs.fam"), quote = F, sep = "\t", row.names = F, col.names = F)
# Create a file for updating phenotypes:
Pheno_info <- FILTERED.MPBC[,c("IID_new","IID_new", "V6")]
write.table(Pheno_info, file=paste("Pheno_info.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# Update sex (--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing) in the (n+2)th column. If no second parameter is provided, n defaults to 1. It is frequently useful to set n=3, since sex defaults to the 5th column in .ped and .fam files.):

plink --bfile MPBC.Rsq03 --update-sex UpdatedIDs.fam 3 --make-bed --out MPBC.Rsq03_sex

# Update phenotype (--make-pheno. If the file has a third column, and a value other than '*' is given, --make-pheno will designate all samples with third column entry equal to the given value as cases (so 2), all other samples mentioned in the file as controls, and all samples missing from the file as having missing phenotypes.):

plink --bfile MPBC.Rsq03_sex --make-pheno Pheno_info.txt 2 --make-bed --out MPBC.Rsq03_sex_pheno

# Filter MAF>5%:
plink --bfile MPBC.Rsq03_sex_pheno --maf 0.05 --make-bed --out MPBC.Rsq03_sex_pheno_MAF05
