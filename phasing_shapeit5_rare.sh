module load StdEnv/2020 r/4.1.2 plink/2.00a3.6 gcc/9.3.0
module load htslib/1.10.2
module load vcftools/0.1.16 bcftools/1.16

# Software directories
shapeit5_common=
shapeit5_rare=
duohmm=
# Pathe and variables
data_dir=
gen_map=
mask=without_mask
chr=$1
#The chunks that we use are provided by the authors of ShapeIt5. It is based on the UKBB data; These chunks only help to parallelize the phasing of rare variants
chunks=/shapeit5/resources/chunks/b38/UKB_WGS_200k/
path_chunks=${data_dir}chunks/
seq_path_data=seq_FINAL_chr${chr}_${mask}.vcf.gz

### Generate bfiles
plink2 --vcf ${seq_path_data} --chr ${chr} --make-bed --keep-allele-order --out ${data_dir}seq_FINAL_chr${chr}_${mask}
### Add genetic positions. Be sure that every markers are in the genetic map, if not, trimming is needed to avoid errors
Rscript gen_pos_interpolation_grch38.R -bim ${data_dir}seq_FINAL_chr${chr}_${mask}.bim -ref_dir ${gen_map} -method "linear"
### Rearrange genetic maps for ShapeIt5
echo -e "pos\tchr\tcM" > ./genmap/plink.chr${chr}.GRCh38.map
awk 'BEGIN{OFS="\t"} {print $4, $1, $3}' ${gen_map}/plink.chr${chr}.GRCh38.map >> ./genmap/plink.chr${chr}.GRCh38.map
### Remove variants not in the map
cat ${gen_map}/plink.chr${chr}.GRCh38.map | cut -d " " -f 1,4 | sed -n '1p;$p' > ${data_dir}genmap/gen_pos_ranges_chr${chr}_tmp.txt
awk '{a[$1] = a[$1] " " $2} END {for (i in a) print i a[i]}' ${data_dir}genmap/gen_pos_ranges_chr${chr}_tmp.txt | awk 'BEGIN {OFS="\t"} {print $0, $1}' > ${data_dir}genmap/gen_pos_ranges_chr${chr}.txt;
rm ${data_dir}genmap/gen_pos_ranges_chr${chr}_tmp.txt
plink2 --bfile ${data_dir}seq_FINAL_chr${chr}_${mask} --extract range ${data_dir}genmap/gen_pos_ranges_chr${chr}.txt --make-bed --keep-allele-order -out ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}
### Generate the .vcf versions without the variants
plink2 --bfile ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr} --recode vcf --keep-allele-order -out ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}

### Match sample IDs in the vcf and fam files (There was a mismatch between the files, IDs repeated during conversion to vcf) 
bcftools query -l ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}.vcf > ${data_dir}samples_chr${chr}.txt
Rscript modify_sample_IDs.R -list ${data_dir}samples_chr${chr}.txt -list_out ${data_dir}new_samples_chr${chr}.txt -fam ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}.fam -fam_out ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_names.fam
bcftools reheader -s ${data_dir}new_samples_chr${chr}.txt -o ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}.vcf

### Phase common variants
echo "Phase common variants"
bcftools +fill-tags ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf | bgzip -c > ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf.gz
bcftools index -t -f ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf.gz
### Reformat genetic map
echo -e "pos chr cM" > ${data_dir}genmap/plink.chr${chr}.GRCh38.fixed.map
awk '{print $4, $1, $3}' ${gen_map}/plink.chr${chr}.GRCh38.map >> ${data_dir}genmap/plink.chr${chr}.GRCh38.fixed.map
### Start phasing
${shapeit5_common} --input ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf.gz \
                   --filter-maf 0.001 \
                   --pedigree ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_names.fam \
                   --region ${chr} \
                   --map ${data_dir}genmap/plink.chr${chr}.GRCh38.fixed.map \
                   --output ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_scaffold.bcf \
                   --thread 8

### Run the phasing of rare variants by chunks, as proposed by the authors
echo "Phase rare variants"
cat ${chunks}/large_chunks_25cM/chunks_chr${chr}.txt > ${path_chunks}/tmp_chunks_chr${chr}.txt
sed "s/chr//g" ${path_chunks}/tmp_chunks_chr${chr}.txt > ${path_chunks}/chunks_chr${chr}.txt
rm ${path_chunks}/tmp_chunks_chr${chr}.txt
while read LINE; do
       CHK=$(echo $LINE | awk '{ print $1; }')
       SRG=$(echo $LINE | awk '{ print $3; }')
       IRG=$(echo $LINE | awk '{ print $4; }')
       ${shapeit5_rare} --input ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_name.vcf.gz \
                   --scaffold ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_scaffold.bcf \
                   --pedigree ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_updated_names.fam \
                   --map ${data_dir}genmap/plink.chr${chr}.GRCh38.fixed.map \
                   --input-region ${IRG} \
                   --scaffold-region ${SRG} \
                   --output ${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_chunk${CHK}.bcf \
                   --thread 8
done < ${path_chunks}/chunks_chr${chr}.txt                   

### Merge the phased files
cut -f 1 ${chunks}/large_chunks_25cM/chunks_chr${chr}.txt | sort -n > ${data_dir}ordered_chunks_ids_chr${chr}.txt
> ${data_dir}chunks.files_chr${chr}.txt
while read CHK; do
  echo "${data_dir}seq_FINAL_${mask}_in_genmap_chr${chr}_chunk${CHK}.bcf" >> ${data_dir}chunks.files_chr${chr}.txt
done < ${data_dir}ordered_chunks_ids_chr${chr}.txt

bcftools concat -n -Ob -o ${data_dir}seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf -f ${data_dir}chunks.files_chr${chr}.txt
bcftools index ${data_dir}seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf
