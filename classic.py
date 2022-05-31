#!/usr/bin/env python
# coding: utf-8


import subprocess as sp
import pandas as pd
import numpy as np
import sys,os,argparse,tempfile,shutil
from os.path import isfile, join


parser = argparse.ArgumentParser(description='Given the fastq files, this program do alignmnt, variant calling and filteration and classify variants with ML algorithm')
parser.add_argument('-r1', '--fastq1', type=str, metavar='<file>', default=None, help='Fastq file first pair')
parser.add_argument('-r2', '--fastq2', type=str, metavar='<file>', default=None, help='Fastq file second pair')
parser.add_argument('-p', '--tools_path', type=str, metavar='<path>', default=None,
                    help='Path to the tools directory', required=True)
parser.add_argument('-w', '--working_dir', type=str, metavar='<path>', default=None,
                    help='the path to the output', required=True)
parser.add_argument('-o', '--output', type=str, default=None, help='Name of the file for writing result')
parser.add_argument('-m', '--running_mode', type=str, default=None, help='Running mode WES or WGS')
parser.add_argument('-t', '--threads', type=str, default=None, help='Number of threads CPU')
parser.add_argument('-r', '--refere.genome_analyzer import GenomeAnalyzernce_file', type=str, metavar='<file>',
                                   help='path to a FASTA-formatted reference file and indexes')

args = parser.parse_args()
tools_path = args.tools_path


if not os.path.exists(args.working_dir + "/QC"):
    os.mkdir(args.working_dir + "/QC")

if not os.path.exists(args.working_dir + "/BAM"):
    os.mkdir(args.working_dir + "/BAM")
    
if not os.path.exists(args.working_dir + "/VCF"):
    os.mkdir(args.working_dir + "/VCF")

if not os.path.exists(args.working_dir + "/Annotation"):
    os.mkdir(args.working_dir + "/Annotation")
    
if not os.path.exists(args.working_dir + "/Final_variants"):
    os.mkdir(args.working_dir + "/Final_variants")

welcomeMessage = """=========================================================
SinglePy.
Analyze WES and WGS raw data
v. 1.0
SinglePy is free non-commercial software. 
Users need to obtain the ANNOVAR licence by themselves. 
Contact the Authors for commercial use.
=============================================================================
"""

endMessage = """============================
Thanks for using SinglePy!
============================================
"""


print (welcomeMessage)



# Alignment and Sorting
fastq_1 = args.fastq1
fastq_2 = args.fastq2
reference = args.reference_file

if None not in (fastq_1, fastq_2, reference):
        bwa_command = "bwa mem" + " " + "-t" + " " + args.threads + " " + args.reference_file + " " + fastq_1 + " " + fastq_2 + " " + | + " " + "java -Xmx10g -jar" + "tools/picard.jar SortSam VALIDATION_STRINGENCY=SILENT" + " " + "O=" + args.working_dir + "BAM/" + args.output + "_hg19_sorted.bam"
        print ("Launching BWA!\n")
else:
        print (args.output + " is not an expected fastq file... skipped...")
        continue

# Running BWA!
    
print ("BWA Command: " + bwa_command + "\n")
process = sp.Popen(bwa_command, shell=True)
process.wait()
print ("Alignment DONE!\n")



# Running Markduplicates
print ("Removing duplicate reads... ")
bwa_out = args.working_dir + "BAM/" + args.output + "_hg19_sorted.bam"

if None not in (bwa_out):
    mkdup_command = "java -jar" + " " + tools_path +"picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT AS=true REMOVE_DUPLICATES=true " + "I=" + bwa_out + " " + "O=" + args.working_dir + "BAM/" + args.output + "_hg19_sorted_mkdup.bam" + " " + "M=" + args.working_dir + "BAM/" + args.output + ".metrics" 
    process = sp.Popen(mkdup_command, shell=True)
    process.wait()
    print ("Duplicates Removed!\n")
else:
    print (args.output + "_hg19_sorted.bam" + " was not found... skipped...")
    continue

    

# FixReadGroup
print ("Fixing  Readgroups... ")
mkdup_out = args.working_dir + "BAM/" + args.output + "_hg19_sorted_mkdup.bam"

if None not in (mkdup_out):
    Fix_command = "java -jar" + " " + tools_path + "picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT" + "I=" + mkdup_out + " " + "O=" + args.working_dir + "BAM/" + args.output + "_hg19_rg.bam" + " " + "RGID=SRR" + atgs.output + " " + "RGLB=SRR" + args.output + " " + "RGPL=illumina" + " " + "RGPU=SRR" + args.output + "RGSM=SRR" + args.output + " " + "CREATE_INDEX=true" + " "+ "&&" + " " + "samtools index" + " " + args.working_dir + "BAM/" + args.output + "_hg19_rg.bam" 
    "CREATE_INDEX=true"
    process = sp.Popen(Fix_command, shell=True)
    process.wait()
    print ("Duplicates Removed!\n")
else:
    print (args.output + "_hg19_sorted_mkdup.bam" + " was not found... skipped...")
    continue


# BAM QC commands


# VariantCalling GATK (haplotypecaller)
print ("VariantCalling with GATK... ")
final_bam = args.working_dir + "BAM/" + args.output + "_hg19_rg.bam"
final_bam_index = args.working_dir + "BAM/" + args.output + "_hg19_rg.bam.bai"

if None not in (final_bam, final_bam_index):
    GATK_command = "java -jar" + " " + tools_path +"gatk-package-4.2.4.0-local.jar HaplotypeCaller" + "-R" + " " + args.reference_file + " " + "-I" + " " + final_bam + "-O" + " " + args.working_dir + "VCF/" + args.output + "_GATK.vcf.gz"
    process = sp.Popen(GATK_command, shell=True)
    process.wait()
    print ("GATK variantcalling completed!\n")
else:
    print (args.output + "_hg19_rg.bam" + " was not found... skipped...")
    continue


# Select SNPs and INDELs
print ("Select SNPs and INDELs... ")
GATK_out = args.working_dir + "VCF/" + args.output + "_GATK.vcf.gz"

if None not in (final_bam, final_bam_index):
    SNP_command = "vcftools --gzvcf" + " " + GATK_out + " " + " --remove-indels --recode --recode-INFO-all" + " " + "--out" + " " + args.working_dir + "VCF/" + args.output + "_GATK_SNP" + " " +  "&&" + " " +  "mv" + " " + args.working_dir + "VCF/" + args.output + "_GATK_SNP.recode.vcf" + " " + args.working_dir + "VCF/" + args.output + "_GATK_SNP.vcf"
    process = sp.Popen(SNP_command, shell=True)
    process.wait()
    indel_command = "vcftools --gzvcf" + " " + GATK_out + " " + " --keep-only-indels --recode --recode-INFO-all" + " " + "--out" + " " + args.working_dir + "VCF/" + args.output + "_GATK_indel" + " " +  "&&" + " " +  "mv" + " " + args.working_dir + "VCF/" + args.output + "_GATK_indel.recode.vcf" + " " + args.working_dir + "VCF/" + args.output + "_GATK_indel.vcf"
    process = sp.Popen(indel_command, shell=True)
    process.wait()
    print ("SNPs & INDELs selected!\n")
else:
    print (args.output + "_hg19_rg.bam" + " was not found... skipped...")
    continue


# Variant hard filtering 
print ("Variant hard filtering... ")
SNP_vcf = args.working_dir + "VCF/" + args.output + "_GATK_SNP.vcf"
indel_vcf = args.working_dir + "VCF/" + args.output + "_GATK_indel.vcf"

if None not in (SNP_vcf, indel_vcf):
    filter_snp_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar VariantFiltration" + " " + "-R" + " " + args.reference_file + "-V" + " " + SNP_vcf + " " + "-O" + " " + args.working_dir + "VCF/" + args.output + "_GATK_SNP_FP.vcf" + " " + '--filter-expression "QD < 2.00" 		--filter-name "QDlessthan2" --filter-expression "FS > 60.000" 		--filter-name "FSgreaterthan60" --filter-expression "SOR > 3.000" 		--filter-name "SORgreaterthan3" --filter-expression "MQRankSum<-2.50" 		--filter-name "MQRankSumlessthanminus2.5" --filter-expression "MQRankSum > 2.50" 		--filter-name "MQRankSumgreaaterthan2.5" --filter-expression "ReadPosRankSum < -1.00" 		--filter-name "ReadPosRankSumlessthanminus1" --filter-expression "ReadPosRankSum > 3.50" 		--filter-name "ReadPosRankSumgreaterthan3.5" --filter-expression "MQ < 60.00" 		--filter-name "MQlessthan60"'
    process.wait()
    filter_indel_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar VariantFiltration" + " " + "-R" + " " + args.reference_file + "-V" + " " + indel_vcf + " " + "-O" + " " + args.working_dir + "VCF/" + args.output + "_GATK_indel_FP.vcf" + " " + '--filter-expression "QD < 2.00" --filter-name "QDlessthan2" 		--filter-expression "FS > 60.000" --filter-name "FSgreaterthan60" 		--filter-expression "SOR > 3.000" --filter-name "SORgreaterthan3" 		--filter-expression "ReadPosRankSum < -1.00" --filter-name "ReadPosRankSumlessthanminus1" 		--filter-expression "ReadPosRankSum > 3.5" --filter-name "ReadPosRankSumgreaterthan3.5" 		--filter-expression "MQ < 60.00" --filter-name "MQlessthan60"'
    process.wait()
    print ("SNPs & INDELs selected!\n")
else:
    print (args.output + "SNP or indel vcf files" + " were not found... skipped...")
    continue


# Variant Recalibration Step
print ("Variant Recalibration Step... ")
SNP_FP = args.working_dir + "VCF/" + args.output + "_GATK_SNP_FP.vcf"
indel_FP = args.working_dir + "VCF/" + args.output + "_GATK_indel_FP.vcf"
miles = 
dbsnp = 
hapmap = 
omni = 
thousand_G = 
if None not in (SNP_FP, indel_FP, miles, dbsnp, omni, thousand_G):
    filter_snp_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar VariantRecalibrator" + " " + "-V" + " " + SNP_FP + " " + '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0     		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP     		-mode SNP     		--max-gaussians 6             -resource:hapmap,known=false,training=true,truth=true,prior=15:' + hapmap + " " + '-resource:omni,known=false,training=true,truth=true,prior=12:' + omni + " " + '-resource:1000G,known=false,training=true,truth=false,prior=10:' + thousand_G + " " + '-resource:dbsnp,known=true,training=false,truth=false,prior=7:' + dbsnp + " " + "-O" + " " + args.working_dir + "VCF/" + args.output + "_snp.recal" + "--tranches-file" + " " + args.working_dir + "VCF/" + args.output +  "_snp.tranches"
    process.wait()
    filter_indel_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar VariantRecalibrator" + " " + "-V" + " " + indel_FP + " " + '--trust-all-polymorphic     	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 		-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 		-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0     	-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \ 
        -mode INDEL     	--max-gaussians 4         -resource:mills,known=false,training=true,truth=true,prior=12:' + miles + " " + "-resource:dbsnp,known=true,training=false,truth=false,prior=2:" + dbsnp + "-O" + " " + args.working_dir + "VCF/" + args.output + "_indels.recal" + " " + "--tranches-file" + " " + args.working_dir + "VCF/" + args.output + "_INDEL.tranches"
    process.wait()
    print ("SNPs & INDELs selected!\n")
else:
    print (args.output + "SNP or indel vcf files and databases" + " were not found... skipped...")
    continue


# ApplyVQSR
print ("Applying VQSR... ")
recalib_SNP = args.working_dir + "VCF/" + args.output + "_snp.recal"
recalib_indel = args.working_dir + "VCF/" + args.output + "_indels.recal"
tranches_snp = args.working_dir + "VCF/" + args.output +  "_snp.tranches"
tranches_indel =  args.working_dir + "VCF/" + args.output + "_INDEL.tranches"
.working_dir + "VCF/" + args.output + "_INDEL.tranches"

if None not in (recalib_SNP, recalib_indel, tranches_snp, tranches_indel):
    VQSR_snp_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar ApplyVQSR" + " " + "-V" + " " + SNP_FP + " " + "--recal-file" + " " + recalib_SNP + " " +  "--tranches-file" + " " + tranches_snp  + " " + "--truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP" + " " + "-O" + " " + args.working_dir + "VCF/" + args.output + "_SNP_recalibrated.vcf.gz"
    process = sp.Popen( VQSR_snp_command, shell=True)
    process.wait()
    VQSR_indel_command = "java -jar" + " " + tools_path + "gatk-package-4.2.4.0-local.jar ApplyVQSR" + " " + "-V" + " " + indel_FP + " " +  "--recal-file" + " " + recalib_indel + " " +  "--tranches-file" + " " + tranches_indel  + " " + "--truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode INDEL" + " " +  "-O" + args.working_dir + "VCF/" + args.output + "_indel_recalibrated.vcf.gz"
    process = sp.Popen(VQSR_indel_command, shell=True)
    process.wait()
    print ("ApplyVQSR Done!\n")
else:
    print (" Recalib and tranche files were not found... skipped...")
    continue

# Call DeepVariant 
print ("VariantCalling with DeepVariant... ")
final_bam = args.working_dir + "BAM/" + args.output + "_hg19_rg.bam"
final_bam_index = args.working_dir + "BAM/" + args.output + "_hg19_rg.bam.bai"

if None not in (final_bam, final_bam_index):
    Deep_command = "singularity run -B /usr/lib/locale/:/usr/lib/locale/" + " " + "docker://google/deepvariant:1.3.1" + " " + "/opt/deepvariant/bin/run_deepvariant" + " " + "--model_type=" + args.running_mode + " " + "--ref=" + args.reference_file + " " + "--reads=" + final_bam  + " " + "--output_vcf=" + args.working_dir + "VCF/" + args.output + "_deepVariant.vcf.gz" + " " + "--intermediate_results_dir" + " " + args.working_dir + "VCF/" + " " + "--num_shards=" + args.threads
    process = sp.Popen(Deep_command, shell=True)
    process.wait()
    print ("DeepVariant Call Done!\n")
else:
    print (" Bam and its index was not found... skipped...")
    continue


# Select SNPs & indels for deepvariant
print ("DeepVariant ouput vcf processing... ")
Deep_vcf = args.working_dir + "VCF/" + args.output + "_deepVariant.vcf.gz"

if None not in (Deep_vcf):
    Deep_SNP_command = "vcftools --gzvcf" + " " + Deep_vcf + " " + " --remove-indels --recode --recode-INFO-all" + " " + "--out" + " " + args.working_dir + "VCF/" + args.output + "_deep_SNP" + " " +  "&&" + " " +  "mv" + " " + args.working_dir + "VCF/" + args.output + "_deep_SNP.recode.vcf" + " " + args.working_dir + "VCF/" + args.output + "_deep_SNP.vcf"
    process = sp.Popen(Deep_SNP_command, shell=True)
    process.wait()
    Deep_indel_command = "vcftools --gzvcf" + " " + Deep_vcf + " " + " --keep-only-indels --recode --recode-INFO-all" + " " + "--out" + " " + args.working_dir + "VCF/" + args.output + "_deep_indel" + " " +  "&&" + " " +  "mv" + " " + args.working_dir + "VCF/" + args.output + "_deep_indel.recode.vcf" + " " + args.working_dir + "VCF/" + args.output + "_deep_indel.vcf"
    process = sp.Popen(Deep_indel_command, shell=True)
    process.wait()
    print ("DeepVariant SNPs & INDELs selected!\n")
else:
    print (args.output + "_deepVariant.vcf.gz" + " was not found... skipped...")
    continue


# Merging vcf files 
print ("Merging vcf files...")
GATK_snp = args.working_dir + "VCF/" + args.output + "_SNP_recalibrated.vcf.gz"
GATK_indel = args.working_dir + "VCF/" + args.output + "_indel_recalibrated.vcf.gz"
Deep_snp = args.working_dir + "VCF/" + args.output + "_deep_SNP.vcf"
Deep_indel = args.working_dir + "VCF/" + args.output + "_deep_indel.vcf"

if None not in (GATK_snp, GATK_indel, Deep_snp, Deep_indel):
    Merge_snp_command = "bcftools merge" + " " + GATK_snp + " " +  Deep_snp + " " + '>' + " " + args.working_dir + "VCF/" + args.output + "_merged_snp.vcf"
    process = sp.Popen(Deep_SNP_command, shell=True)
    process.wait()
    Merge_indel_command = "bcftools merge" + " " + GATK_indel + " " + Deep_indel + " " + '>' + " " + args.working_dir + "VCF/" + args.output + "_merged_indel.vcf"
    process = sp.Popen(Deep_indel_command, shell=True)
    process.wait()
    print ("Vcf files merged!\n")
else:
    print ("GATK, DeepVariant vcf files" + " were not found... skipped...")
    continue

