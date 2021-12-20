#!/bin/bash

module load bcftools
module load tabix/default


##INPUTS
vcf_in=$1
vcf_name=$(basename $vcf_in | sed 's/.gz//' | sed 's/.vcf//')
out_dir=$2
toml="/data/talkowski/ev962/cffDNA/cffDNA_analysis/annotation/cffDNA_genelist_missense_config.toml"

#### Annotate VCF with annovar ####
perl /PHShome/ev962/software/annovar_20200608/table_annovar.pl  $vcf_in  /PHShome/ev962/software/annovar_20200608/humandb/ -otherinfo -buildver hg38 -out $out_dir/$vcf_name"_annovar" -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad211_exome,gnomad30_genome,avsnp150,dbnsfp41a,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,clinvar_20210501,gene4denovo201907 -operation g,g,r,f,f,f,f,f,r,f,f,f,f -nastring . -vcfinput

bgzip $out_dir/$vcf_name"_annovar.hg38_multianno.vcf"

#### Annotate with gene list ###
##Use vcfanno to annotate with disease gene list and regional missense constraint. If want to change lists, edit the *toml file
/PHShome/ev962/software/vcfanno_linux64 -base-path ${out_dir} ${toml} ${out_dir}/${vcf_name}"_annovar.hg38_multianno.vcf.gz" | bgzip > ${out_dir}/${vcf_name}"_annovar.hg38_multianno_genelist.vcf.gz"

tabix -p vcf $${out_dir}/${vcf_name}"_annovar.hg38_multianno_genelist.vcf.gz"

#### Remove old vcfs and files needed only for annotation to save space ####
rm ${out_dir}/${vcf_name}"_annovar.hg38_multianno.vcf.gz"
rm ${out_dir}/${vcf_name}"_annovar.avinput"
rm ${out_dir}/${vcf_name}"_annovar.hg38_multianno.txt"
