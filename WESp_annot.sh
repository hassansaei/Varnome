
r1=M026

# Directories
INPUT_DIR="${PWD}/quickstart-input"
OUTPUT_DIR="${PWD}/quickstart-output"
INPUT_REF="${PWD}/quickstart-ref"

##############################
# Adding annotation with ANNOVAR
##############################
OUTPUT_DIR="${PWD}/quickstart-output"
cd quickstart-output/

perl  ../bin/InterVar-master/convert2annovar.pl -format vcf4 "$r1"_Final_SNP.vcf -outfile "$r1"_Final_SNP.avinput 
perl  ../bin/InterVar-master/convert2annovar.pl -format vcf4 "$r1"_Final_INDEL.vcf -outfile "$r1"_Final_INDEL.avinput 

perl ../bin/InterVar-master/table_annovar.pl "$r1"_Final_SNP.avinput ../bin/InterVar-master/humandb/ -buildver hg19 -out "$r1"_SNP -remove -protocol refGene,esp6500siv2_all,avsnp150,dbnsfp41a,clinvar_20210123,gnomad211_genome,dbscsnv11,dbnsfp31a_interpro,popfreq_all_20150413,exac03,ensGene,knownGene -operation g,f,f,f,f,f,f,f,f,f,g,g -nastring . -polish -xref gene_xref.txt --otherinfo

perl ../bin/InterVar-master/table_annovar.pl "$r1"_Final_INDEL.avinput ../bin/InterVar-master/humandb/ -buildver hg19 -out "$r1"_INDEL -remove -protocol refGene,esp6500siv2_all,avsnp150,dbnsfp41a,clinvar_20210123,gnomad211_genome,dbscsnv11,dbnsfp31a_interpro,popfreq_all_20150413,exac03,ensGene,knownGene -operation g,f,f,f,f,f,f,f,f,f,g,g -nastring . -polish -xref gene_xref.txt --otherinfo


#######################################
# ACMG-AMP interpertation with InterVar
#######################################
OUTPUT_DIR="${PWD}/quickstart-output"

python3 ../bin/InterVar-master/Intervar.py -b hg19 -i "$r1"_SNP.hg19_multianno.txt  -d ../bin/InterVar-master/humandb/ -t ../bin/InterVar-master/intervardb/  --skip_annovar -o "$r1"_SNP

python3 ../bin/InterVar-master/Intervar.py -b hg19 -i "$r1"_INDEL.hg19_multianno.txt  -d ../bin/InterVar-master/humandb/ -t ../bin/InterVar-master/intervardb/  --skip_annovar -o "$r1"_INDEL
