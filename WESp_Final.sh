
read -p " Please enter your reference:" Re1
read -p " Please enter the read 1:" r1
read -p " Please enter the read 2:" r2
read -p " Number of Threads:" T
# Directories

INPUT_DIR="${PWD}/quickstart-input"
OUTPUT_DIR="${PWD}/quickstart-output"
INPUT_QC="${PWD}/bin/FastQC"
INPUT_REF="${PWD}/quickstart-ref/BWA_index"

########
# FASTQC
########
echo "Data Quality check"

"${INPUT_QC}"/./fastqc "${INPUT_DIR}"/*.fq.gz

###################
# Mapping/alignment
###################

echo "Mapping/alignment with BWA"
cd "${INPUT_DIR}"
bwa mem -t $T -R "@RG\tID:A8100\tSM:A8100" "${INPUT_REF}"/"$Re1".fa "${INPUT_DIR}"/"$r1"_1.fq.gz "${INPUT_DIR}"/"$r2"_2.fq.gz -o "$r1".sam 


#############
# SAM to BAM
############
echo "SAM to BAM and Sort"
mv *.sam ../quickstart-output/
cd ../quickstart-output/

samtools view -@ $T -bS "$r1".sam > "$r1".bam
bamtools stats -insert -in "$r1".bam > "$r1"-summary
java -Xmx10g -jar ../bin/picard.jar SortSam VALIDATION_STRINGENCY=SILENT I="$r1".bam O="$r1"_Sort.bam SORT_ORDER=coordinate
rm *.sam

##########
# Markdup
##########
echo "Remove Duplicates"

java -jar ../bin/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT AS=true REMOVE_DUPLICATES=true I="$r1"_Sort.bam O="$r1"_sort_mkdup.bam M="$r1".metrics

java -jar ~/bin/gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar DepthOfCoverage -R ~/quickstart-ref/BWA_index/hg19_ucsc.fa -O result -I RZL.hg19-bwamem.bam -gene-list All_gene.refseq -L truseq-dna-exome-targeted-regions-manifest-v1-2.bed



########################
# AddorReplaceReadGroup
#######################
echo "ADD ReadGroup"

java -jar ../bin/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I="$r1"_sort_mkdup.bam O="$r1"_sort_mkdup_rg.bam  RGID=SRR"$r1" RGLB=SRR"$r1" RGPL=illumina RGPU=SRR"$r1" RGSM=SRR"$r1" CREATE_INDEX=true

samtools view -H "$r1"_sort_mkdup_rg.bam | grep '^@RG'

###########
# index BAM
###########
echo "BAM file indexing"

samtools index "$r1"_sort_mkdup_rg.bam

rm *_Sort.bam *_sort_mkdup.bam *.metrics "$r1".bam

#gzip  "$r1"_sort_mkdup_rg.bam
#gzip  "$r1"_sort_mkdup_rg.bam.bai


########################
# Variant Calling GATK
#######################
cd quickstart-output/
echo "Variant Calling  with GATK"

java -jar ../bin/gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar HaplotypeCaller -R "${INPUT_REF}"/"$Re1".fa -I "${OUTPUT_DIR}"/"$r1"_sort_mkdup_rg.bam -O "${OUTPUT_DIR}"/"$r1"_GATK.vcf.gz


########################
# Seperate SNV & INDELS
#######################
echo "Seperate SNV & INDELS"

vcftools --gzvcf "$r1"_GATK.vcf.gz --remove-indels --recode --recode-INFO-all --out "$r1"_GATK_SNP
vcftools --gzvcf "$r1"_GATK.vcf.gz --keep-only-indels  --recode --recode-INFO-all --out "$r1"_GATK_INDEL

mv "$r1"_GATK_SNP.recode.vcf ./"$r1"_GATK_SNP.vcf
mv "$r1"_GATK_INDEL.recode.vcf ./"$r1"_GATK_INDEL.vcf

java -jar ../bin/gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar VariantFiltration -R "${INPUT_REF}"/"$Re1".fa -V "$r1"_GATK_SNP.vcf -O "$r1"_GATK_SNP_F.vcf --filter-expression "QD < 2.00" --filter-name "QDlessthan2" --filter-expression "FS > 60.000" --filter-name "FSgreaterthan60" --filter-expression "SOR > 3.000" --filter-name "SORgreaterthan3" --filter-expression "MQRankSum<-2.50" --filter-name "MQRankSumlessthanminus2.5" --filter-expression "MQRankSum > 2.50" --filter-name "MQRankSumgreaaterthan2.5" --filter-expression "ReadPosRankSum < -1.00" --filter-name "ReadPosRankSumlessthanminus1" --filter-expression "ReadPosRankSum > 3.50" --filter-name "ReadPosRankSumgreaterthan3.5" --filter-expression "MQ < 60.00" --filter-name "MQlessthan60"

java -jar ../bin/gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar VariantFiltration -R "${INPUT_REF}"/"$Re1".fa -V "$r1"_GATK_INDEL.vcf -O "$r1"_GATK_INDEL_F.vcf --filter-expression "QD < 2.00" --filter-name "QDlessthan2" --filter-expression "FS > 60.000" --filter-name "FSgreaterthan60" --filter-expression "SOR > 3.000" --filter-name "SORgreaterthan3" --filter-expression "ReadPosRankSum < -1.00" --filter-name "ReadPosRankSumlessthanminus1" --filter-expression "ReadPosRankSum > 3.5" --filter-name "ReadPosRankSumgreaterthan3.5" --filter-expression "MQ < 60.00" --filter-name "MQlessthan60"

grep -E '^#|PASS' "$r1"_GATK_SNP_F.vcf > "$r1"_GATK_SNP_FP.vcf
grep -E '^#|PASS' "$r1"_GATK_INDEL_F.vcf > "$r1"_GATK_INDEL_FP.vcf

grep -vc "^#" "$r1"_GATK_SNP_F.vcf 
grep -vc "^#" "$r1"_GATK_SNP_FP.vcf 

grep -vc "^#" "$r1"_GATK_INDEL_F.vcf 
grep -vc "^#" "$r1"_GATK_INDEL_FP.vcf 

#####################################################
# Variant Calling  with AI and SNV & INDEL separation
####################################################

cd quickstart-output/
echo "Variant Calling  with AI step 1"
N_SHARDS=10
LOGDIR=./BWA_logs
mkdir -p "${LOGDIR}"

time seq 0 $((N_SHARDS-1)) | parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" python ../bin/deepvariant/binaries/DeepVariant-1.0.0/make_examples.zip  --mode calling --ref "${INPUT_REF}"/"$Re1".fa --reads "$r1"_sort_mkdup_rg.bam --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" --task {}

python ../bin/deepvariant/binaries/DeepVariant-1.0.0/call_variants.zip --outfile "$r1"_Deep.gz --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" --checkpoint ../bin/deepvariant/models/DeepVariant/1.0.0/DeepVariant-inception_v3-1.0.0+data-wgs_standard/model.ckpt

python ../bin/deepvariant/binaries/DeepVariant-1.0.0/postprocess_variants.zip --ref "${INPUT_REF}"/"$Re1".fa --infile "$r1"_Deep.gz --outfile "$r1"_Deep.vcf


vcftools --vcf "$r1"_Deep.vcf --remove-indels --recode --recode-INFO-all --out "$r1"_Deep_SNP
vcftools --vcf "$r1"_Deep.vcf --keep-only-indels  --recode --recode-INFO-all --out "$r1"_Deep_INDEL

mv "$r1"_Deep_SNP.recode.vcf ./"$r1"_Deep_SNP.vcf
mv "$r1"_Deep_INDEL.recode.vcf ./"$r1"_Deep_INDEL.vcf

rm examples.*.gz
rm *.log
rm *.sam

'''
##############################
# Adding annotation with ANNOVAR
##############################
OUTPUT_DIR="${PWD}/quickstart-output"
cd quickstart-output/

perl  ../bin/InterVar-master/convert2annovar.pl -format vcf4 "$r1"_Final_SNP.vcf -outfile "$r1"_Final_SNP.avinput 
perl  ../bin/InterVar-master/convert2annovar.pl -format vcf4 "$r1"_Final_INDEL.vcf -outfile "$r1"_Final_INDEL.avinput 

perl ../bin/Annovar/table_annovar.pl MohammadHossein_Alishahi_GATK_SNP.avinput ../bin/Annovar/humandb/ -buildver hg19 -out MohammadHossein_Alishahi_SNP -remove -protocol refGene,knownGene,ensGene,avsift,dbnsfp42a,intervar_20180118,esp6500siv2_all,exac03,gnomad211_genome,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20210501,cadd13gt10,popfreq_all_20150413   -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -polish -xref gene_xref.txt --otherinfo


perl ../bin/Annovar/table_annovar.pl MohammadHossein_Alishahi_GATK_SNP.avinput ../bin/Annovar/humandb/ -buildver hg19 -out MohammadHossein_Alishahi_SNP -remove -protocol refGene,knownGene,ensGene,avsift,dbnsfp42a,esp6500siv2_all,exac03,gnomad211_genome,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20210501,cadd13gt10,popfreq_all_20150413   -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -polish -xref gene_xref.txt --otherinfo

#######################################
# ACMG-AMP interpertation with InterVar
#######################################
OUTPUT_DIR="${PWD}/quickstart-output"

python3 ../bin/InterVar-master/Intervar.py -b hg19 -i "$r1"_SNP.hg19_multianno.txt  -d ../bin/InterVar-master/humandb/ -t ../bin/InterVar-master/intervardb/  --skip_annovar -o "$r1"_SNP

python3 ../bin/InterVar-master/Intervar.py -b hg19 -i "$r1"_INDEL.hg19_multianno.txt  -d ../bin/InterVar-master/humandb/ -t ../bin/InterVar-master/intervardb/  --skip_annovar -o "$r1"_INDEL

'''

