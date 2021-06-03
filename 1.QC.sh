### Genotyping quality control (QC) in the MPBC GWAS 
### Created by: Kajsa Brolin
### Original code retrieved from: GP2 Training and Networking (https://github.com/GP2-TNC-WG/GP2-Bioinformatics-course). (If you haven't checked it out then please do, it's great!)


# NOTE: do all samples have a gender and a phenotype?

# Make sure that the following programs are installed:
# Bash shell
# PLINK v1.9
# R
# BCFTools
# Perl

# Set your working directory:
cd Desktop/MPBC

### 1) SAMPLE QC
##Calculate genotyping call rates per sample
plink --bfile MPBC --missing --out call_rates
#Remove call rate outliers
plink --bfile MPBC --mind 0.05 --make-bed --out MPBC_after_call_rate

###HETEROZYGOSITY
plink --bfile MPBC_after_call_rate --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile MPBC_after_call_rate --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet

awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
mv prunedHet.het HETEROZYGOSITY_DATA.txt

plink --bfile MPBC_after_call_rate --remove all_outliers.txt --make-bed --out MPBC_after_heterozyg

###GENDER CHECKING
plink --bfile MPBC_after_heterozyg --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
plink --bfile MPBC_after_heterozyg --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 

###ANCESTRY
plink --bfile MPBC_after_heterozyg --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile MPBC_after_heterozyg --flip hapmap3_bin_snplis-merge.missnp --make-bed --out MPBC_after_heterozyg3
plink --bfile MPBC_after_heterozyg3 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile MPBC_after_heterozyg3 --exclude hapmap3_bin_snplis-merge.missnp --out MPBC_after_heterozyg4 --make-bed
plink --bfile MPBC_after_heterozyg4 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile hapmap3_bin_snplis --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile hapmap3_bin_snplis --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet
plink --bfile pruned_data --geno 0.05 –exclude range exclusion_regions_hg19.txt --pca --out pca --make-bed

grep "EUROPE" pca.eigenvec > eur.txt
grep "ASIA" pca.eigenvec > asia.txt
grep "AFRICA" pca.eigenvec > afri.txt
grep -v -f eur.txt pca.eigenvec | grep -v -f asia.txt | grep -v -f afri.txt > new_samples.txt
cut -d " " -f 3 after_gender.fam > new_samples_add.txt
paste new_samples_add.txt new_samples.txt > new_samples2.txt
paste eur_add.txt eur.txt > euro.txt
paste asia_add.txt asia.txt > asiao.txt
paste afri_add.txt afri.txt > afrio.txt
cat new_samples2.txt euro.txt asiao.txt afrio.txt > pca.eigenvec2

# Run 1checkPCA.R to plot and identify ancestry outliers
# Continue with the following after:

plink --bfile MPBC_after_heterozyg --keep europeanAncestrySampleIds.txt --make-bed --out MPBC_Euro
cat asianAncestrySampleIds.txt africanAncestrySampleIds.txt > hapmap_outliers.txt

###RELATEDNESS
plink --bfile MPBC_Euro --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile MPBC_Euro --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet

gcta64 --bfile pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
plink --bfile MPBC_Euro --keep GRM_matrix_0125.grm.id --make-bed --out MPBC_Euro_pihat

cut -f 1,2 MPBC_Euro.fam > IDs_before_relatedness_filter.txt
cut -f 1,2 MPBC_Euro_pihat.fam > IDs_after_relatedness_filter.txt



### 2) VARIANT QC

###MISSINGNESS PER VARIANT
plink --bfile MPBC_Euro_pihat --make-bed --out MPBC_Euro_pihat_mind --geno 0.05
grep "(--geno)" MPBC_Euro_pihat_mind.log > MISSINGNESS_SNPS.txt

##MISSINGNESS BY CASE-CONTROL P > 1E-4 
plink --bfile MPBC_Euro_pihat_mind --test-missing --out missing_snps 
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
plink --bfile MPBC_Euro_pihat_mind --exclude missing_snps_1E4.txt --make-bed --out MPBC_Euro_pihat_mind_missing1
sort -u missing_snps_1E4.txt > VARIANT_TEST_MISSING_SNPS.txt

##MISSINGNESS BY HAPLOTYPE
plink --bfile MPBC_Euro_pihat_mind_missing1 --test-mishap --out missing_hap 
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\/g/' missing_haps_1E4.txt > missing_haps_1E4_final.txt 

sort -u missing_haps_1E4_final.txt > HAPLOTYPE_TEST_MISSING_SNPS.txt
plink --bfile MPBC_Euro_pihat_mind_missing1 --exclude missing_haps_1E4_final.txt --make-bed --out MPBC_Euro_pihat_mind_missing12

##HWE
plink --bfile MPBC_Euro_pihat_mind_missing12 --filter-controls --hwe 1E-4 --write-snplist --out HWE_snps
plink --bfile MPBC_Euro_pihat_mind_missing12 --extract HWE_snps.snplist --make-bed --out MPBC_Euro_pihat_mind_missing123

##MAF and filtering for autosomal chromosomes only
plink --bfile FILTERED.MPBC --maf 0.01 --autosomal --make-bed --out MPBC_MPBC_Euro_pihat_mind_missing123_MAF

##Rename files:
mv MPBC_MPBC_Euro_pihat_mind_missing123_MAF.bim FILTERED.MPBC.bim
mv MPBC_MPBC_Euro_pihat_mind_missing123_MAF.bed FILTERED.MPBC.bed
mv MPBC_MPBC_Euro_pihat_mind_missing123_MAF.fam FILTERED.MPBC.fam



### POST QC FORMATTING BEFORE IMPUTATION ###

### CHECK DATA AGAINST IMPUTATION REFERENCE SNP LIST
## Reference panel used here: Haplotype Reference Consortium (HRC) version r1.1 2016, European population
## Imputation is done at the Michigan imputation server (MIS)
# Using William Rayner's tool: https://www.well.ox.ac.uk/~wrayner/tools/
# Reference list we will use it against: http://www.haplotype-reference-consortium.org/site

# William Rayner script requires frequency output as an input
plink --bfile FILTERED.MPBC --freq --out FILTERED.MPBC
perl HRC-1000G-check-bim.pl -b FILTERED.MPBC.bim -f FILTERED.MPBC.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
# Creates Run-plink.sh
sh Run-plink.sh

#Now we should see FILTERED.test-updated-chr for each chromosomes.

###CONVERT TO VCF, SORT AND COMPRESS
#When converting to VCF, we will split this to different chromosomes as required by Michigan Imputation Server (it also makes files smaller and easier to handle)

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink --bfile FILTERED.MPBC-updated-chr$chnum --recode vcf --chr $chnum --out MPBC$chnum
done

#MIS also requires bgzip compression. BCFTools can compress and sort the VCF with one command.
export PATH=$PATH:/Applications/bcftools-1.12

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  bcftools sort MPBC$chnum.vcf -Oz -o pre_impute_MPBC$chnum.vcf.gz
  vcf-sort MPBC$chnum.vcf | bgzip -c > pre_impute_MPBC$chnum.vcf.gz
done

#UPLOAD THE pre_impute_MPBC$chnum.vcf.gz FILES TO MICHIGAN IMPUTATION SERVER FOR IMPUTATION

