module load StdEnv/2020 r/4.1.2 plink/2.00a3.6 gcc/9.3.0
module load htslib/1.10.2
module load vcftools/0.1.16 bcftools/1.16

### Define paths and variables
shapeit5_common=
gen_map=
mask=without_mask
file_prefix=
wgs_data=
chr=$1
#Results will be exported here:
data_dir=

### Remove probematic and absent positions
awk '$5=="D" || $6=="D" || $5=="I" || $6=="I" || $5==0 || $6==0 || $5=="+" || $6=="+" || $5=="-" || $6=="-" {print $2}' ${data_dir}${file_prefix}.${chr}.bim  > ${data_dir}${file_prefix}.${chr}.snps2remove.txt
awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C")  {print $2}' ${data_dir}${file_prefix}.${chr}.bim  >> ${data_dir}${file_prefix}.${chr}.snps2remove.txt
plink2 --bfile ${data_dir}${file_prefix}.${chr} --exclude ${data_dir}${file_prefix}.${chr}.snps2remove.txt --allow-extra-chr --make-bed --out ${data_dir}${file_prefix}_snpsexcluded_chr${chr}

### Add genetic positions. Be sure that every markers are in the genetic map, if not, trimming is needed to avoid errors.
Rscript gen_pos_interpolation_grch38.R -bim ${data_dir}${file_prefix}_snpsexcluded_chr${chr}.bim -ref_dir ${gen_map} -method "linear"
### Rearrange genetic maps for ShapeIt5
echo -e "pos\tchr\tcM" > ${data_dir}genmap/plink.genos.chr${chr}.GRCh38.map
awk 'BEGIN{OFS="\t"} {print $4, $1, $3}' ${gen_map}/plink.chr${chr}.GRCh38.map >> ${data_dir}genmap/plink.genos.chr${chr}.GRCh38.map

### Remove variants that are not included in the genetic map
cat ${gen_map}/plink.chr${chr}.GRCh38.map | cut -d " " -f 1,4 | sed -n '1p;$p' > ${data_dir}genmap/gen_pos_ranges_genos_chr${chr}_tmp.txt
awk '{a[$1] = a[$1] " " $2} END {for (i in a) print i a[i]}' ${data_dir}genmap/gen_pos_ranges_genos_chr${chr}_tmp.txt | awk 'BEGIN {OFS="\t"} {print $0, $1}' > ${data_dir}genmap/gen_pos_ranges_genos_chr${chr}.txt;
plink2 --bfile ${data_dir}${file_prefix}_snpsexcluded_chr${chr} --extract range ${data_dir}genmap/gen_pos_ranges_genos_chr${chr}.txt --make-bed -out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}
plink2 --bfile ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr} --rm-dup force-first --make-bed --out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}

### Keep variants that are present in the WGS + flip stand and switch  allele order if needed
bcftools query -f '%ID\t%CHROM\t%POS\t%ALT\t%REF\n'  ${wgs_data}${mask}_in_genmap_chr${chr}.bcf | awk '{print $2"\t"$1"\t0\t"$3"\t"$4"\t"$5}' > ${data_dir}seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bim

### Change all variants to chr:pos:a1:a2
plink2 --bfile ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr} --set-all-var-ids chr@:#:\$1:\$2 --make-bed -out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}
Rscript 2023-12-05_bim_matching.R -bim1 ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr} -bim2 ${data_dir}seq_FINAL_PHASED_${mask}_in_genmap_chr${chr} -removeindels 1 -out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}
plink2 --bfile ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr} --rm-dup force-first --make-bed --out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}

### Generate the .vcf versions without the variants not in the sequences
module load StdEnv/2020 plink/1.9b_6.21-x86_64
# Remove duplicated ids from flip file
sort ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toflip.txt | uniq > ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toflip_cleaned.txt
# Remove duplicated ids from switch file
sort ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toswitch.txt | uniq > ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toswitch_cleaned.txt

plink --bfile ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr} \
      --exclude ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toremove.txt \
      --reference-allele ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toswitch_cleaned.txt \
      --flip ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_toflip_cleaned.txt \
      --make-bed -out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}

plink --bfile ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr} \
      --update-name ${data_dir}${file_prefix}_snpsexcluded_in_genmap_chr${chr}_torrename.txt \
      --recode vcf --keep-allele-order \
      --out ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}

### Create fam and map files in the format of shapeit5. In the .fam file, the fam needs to be appended to the IDs
awk '{print $1"_"$2, $1"_"$3, $1"_"$4}' ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.fam | awk -v OFS='\t'  -F' ' '{for(i=1;i<=NF;i++) if($i~/_0$/) $i=0}1' > ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}_shapeit5.fam

### Phase common variants
bcftools +fill-tags ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.vcf | bgzip -c > ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.vcf.gz
bcftools index -t -f ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.vcf.gz
${shapeit5_common} --input ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.vcf.gz \
                   --pedigree ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}_shapeit5.fam \
                   --region ${chr} \
                   --map ${data_dir}genmap/plink.genos.chr${chr}.GRCh38.map \
                   --output ${data_dir}${file_prefix}_snpsexcluded_in_genmap_matched_chr${chr}.bcf \
                   --thread 8
