module load StdEnv/2020 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 impute5/1.1.5 r/4.1.2
module load plink/1.9b_6.21-x86_64

# Paths and variables
path_impu=
gen_map=
mask=without_mask
chr=$1
vcf_file1=
### Chips to merge (outputs from imputation_impute5.sh)
plink_file1=
plink_file2=
plink_file3=
plink_file4=
plink_file5=
plink_file6=

bcftools merge --force-samples ${path_impu}/${plink_file1} ${path_impu}/${plink_file2} ${path_impu}/${plink_file3} ${path_impu}/${plink_file4} ${path_impu}/${plink_file5} ${path_impu}/${plink_file6} -Oz -o ${path_impu}/CaG_imputed_allchips_filtered_chr${chr}.vcf.gz
bcftools index -t -f ${path_impu}/CaG_imputed_allchips_filtered_chr${chr}.vcf.gz
### Eliminate the missing data as they are doubled (imputed + non imputed)
plink --vcf ${path_impu}/CaG_imputed_allchips_filtered_chr${chr}.vcf.gz --missing --out ${path_impu}/CaG_imputed_allchips_filtered_chr${chr}
plink --vcf ${path_impu}/CaG_imputed_allchips_filtered_chr${chr}.vcf.gz --geno 0.01 --make-bed --out ${path_impu}/CaG_imputed_allchips_filtered_geno0.01_chr${chr}

### Remove positions from the chips and put back all REF alleles
python remove_chip_position_from_imputed_vcf.py ${path_impu}/${vcf_file1}.filtered
bgzip ${path_impu}/${vcf_file1}.filtered.onlyIMP.vcf
bcftools index ${path_impu}/${vcf_file1}.filtered.onlyIMP.vcf.gz

echo 'add samples'
python add_inds_in_vcf.py ${path_impu}/CaG_imputed_chr${chr}
echo 'bgzip vcf file'
bgzip ${path_impu}/CaG_imputed_chr${chr}.onlyIMP2.vcf
echo 'index vcf file'
bcftools index ${path_impu}/CaG_imputed_chr${chr}.onlyIMP2.vcf.gz
echo 'produce bim files'
bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\n' ${path_impu}/CaG_imputed_chr${chr}.onlyIMP2.vcf.gz | awk '{print $1"\t"$2"\t0\t"$3"\t"$4"\t"$5}' > ${path_impu}/CaG_imputed_chr${chr}.onlyIMP2.bim
echo 'compare bim files'
wgs_real_ref_file=/wgs_qc_freeze1/chr${chr}_FINAL_realREF.txt
Rscript bim_matching_for_realREF.R -bim1 ${path_impu}/CaG_imputed_chr${chr}.onlyIMP2.bim -bim2 ${wgs_real_ref_file} -removeindels 0
plink --vcf ${path_impu}/CaG_imputed_chr${chr}.onlyIMP.vcf.gz --reference-allele ${path_impu}/CaG_imputed_chr${chr}.onlyIMP_toswitch.txt --make-bed -out ${path_impu}/CaG_imputed_chr${chr}.onlyIMP.realREF

