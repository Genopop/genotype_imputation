module load StdEnv/2020 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 impute5/1.1.5 r/4.1.2
module load plink/1.9b_6.21-x86_64

### Define paths and variables
path_phasing=
path_impu=
plink_file=  # Genotype data
chunks=
path_chunks=
mask=without_mask
chr=$1
cd ${path_impu}

### Compute the recombination rate and generate the genetic map in the correct format
mkdir=${path_impu}/gen_map
if [ ! -s "${path_impu}/gen_map" ]
then
Rscript gen_map_phaseit2_format.R -ref_dir ${gen_map} -out ${path_impu}/gen_map
fi

### Imputation done with chunk
chunk_len=25
ch=large
cat ${chunks}/${ch}_chunks_${chunk_len}cM/chunks_chr${chr}.txt > ${path_chunks}/tmp_chunks_${chunk_len}cM_chr${chr}.txt
sed "s/chr//g" ${path_chunks}/tmp_chunks_${chunk_len}cM_chr${chr}.txt > ${path_chunks}/chunks_${chunk_len}cM_chr${chr}.txt

while read LINE; do
          CHK=$(echo $LINE | awk '{ print $1; }')
          SRG=$(echo $LINE | awk '{ print $3; }')
          IRG=$(echo $LINE | awk '{ print $4; }')
          echo ${CHK}
          echo ${SRG}
          echo ${IRG}
  impute5 --h ${path_phasing}phasing_wgs/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf \
          --g ${path_phasing}phasing_genotypes/${plink_file}_snpsexcluded_in_genmap_matched_chr${chr}.bcf \
          --m ${path_impu}gen_map/plink.chr${chr}.GRCh38.map \
          --l ${path_impu}logs/log_${plink_file}_chr${chr} \
          --r ${IRG} \
          --buffer-region ${SRG} \
          --o ${path_impu}${plink_file}_chr${chr}_chunk${CHK}_IMPUTED.vcf.gz \
          --threads 24
  bcftools index -t -f ${path_impu}${plink_file}_chr${chr}_chunk${CHK}_IMPUTED.vcf.gz
done < ${path_chunks}/chunks_${chunk_len}cM_chr${chr}.txt

### Combine chunks for each chip-chr
chunks_file=$(cut -f 1 ${chunks}/${ch}_chunks_${chunk_len}cM/chunks_chr${chr}.txt)
for CHK in ${chunks_file}; do echo "${path_impu}${plink_file}_chr${chr}_chunk${CHK}_IMPUTED.vcf.gz" >> ${path_impu}${plink_file}.chunks.files_chr${chr}.txt; done
bcftools concat -n -Oz -o ${path_impu}${plink_file}_chr${chr}_IMPUTED.vcf.gz -f ${path_impu}${plink_file}.chunks.files_chr${chr}.txt
bcftools index ${path_impu}${plink_file}_chr${chr}_IMPUTED.vcf.gz

### Get the INFO scores of every variants (Ranges from 0 to 1, 1 being a perfect imputation accuracy)
final_data=${path_impu}${plink_file}_chr${chr}_IMPUTED.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${final_data} > ${final_data}_tmp_1.txt
bcftools query -f '%INFO\n' ${final_data} | tr ";" "\n" | grep INFO > ${final_data}_tmp_2.txt
paste ${final_data}_tmp_1.txt ${final_data}_tmp_2.txt > ${final_data}_INFO.txt
rm ${final_data}_tmp_1.txt
rm ${final_data}_tmp_2.txt

### Look at switch and flips; 
# remove those which are not imputed at the same positions (the chip non imputed positions are not in the right orientation)
vcf_file=${path_impu}${plink_file}_chr${chr}_IMPUTED
wgs_path=/phasing_wgs/
wgs=seq_FINAL_PHASED_without_mask_in_genmap_chr${chr}
bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\n' ${wgs_path}${wgs}.bcf | awk '{print $1"\t"$2"\t0\t"$3"\t"$4"\t"$5}' > ${wgs}.bim
bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\n' ${vcf_file}.vcf.gz | awk '{print $1"\t"$2"\t0\t"$3"\t"$4"\t"$5}' > ${vcf_file}.bim
Rscript ${path_scripts}2023-12-05_bim_matching.R -bim1 ${vcf_file} -bim2 ${wgs} -removeindels 0 -out ${vcf_file}_WGS
cat ${vcf_file}_WGS_toflip.txt >> ${vcf_file}_WGS_toswitch.txt
cat ${vcf_file}_WGS_toswitch.txt | wc -l
python filter_vcf.py ${vcf_file} ${vcf_file}_WGS_toswitch.txt
bgzip ${vcf_file}.filtered.vcf
bcftools index ${vcf_file}.filtered.vcf.gz




