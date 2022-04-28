version 1.0.0

## Copyright Imagine Institute, 2022
##
## This WDL pipeline implements data preprocessing acording to the GATK Best Practices
##
## Requirements/expectations:
## -
## 
##
## Output :
## - A clean BAM file and its index, ready for variant detection and CNV analysis
## - Tow VCF files and their indexs obtained from GATK HaplotypeCaller and DeepVariant algorithm
## - Annotation files
##
##
## Software version requirments :
## - BWA 
## - GATK
## - Samtools
## - Picard
## - Python
##
##
##

# Workflow Definition
workflow DataPreprocessing {
	String pipeline_version = "1.0.0"
	input {
	String sample_name	
	File r1fastq
	File r2fastq
	File ref_fasta
	File ref_fasta_amb
	File ref_fasta_sa
	File ref_fasta_bwt
	File ref_fasta_ann
	File ref_fasta_pac
	File ref_dict
	File ref_index
	File dbSNO_vcf
	File dbSNP_vcf_index
	File All_gene
	File hapmap_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf
        File omni_resource_vcf_index
        File thousand_G_resource_vcf
        File thousand_G_resource_vcf_index
        File mills_resource_vcf
        File mills_resource_vcf_index
        File dbsnp_vcf
        File dbsnp_vcf_index
	
	Array[File] known_indels_vcf
	Array[File] known_indels_indices
	Array[File] All_genes_hg19
	Array[File] Target_bed 
	Boolean make_gvcf = true
	String gatk_path = "~/Singlypy/tools/gatk/"
	String picard_path = "~/Singlypy/tools/"
	String ref_path_hg19 = "~/Singlepy/References/hg19"
	String ref_path_hg38 = "~/Singlepy/References/hg38"
	String DeepVariant_path = "~/Singlepy/tools/"
	String QC_path = "~/Singlepy/Results/QC/"
	String BAM_path = "~/Singlepy/Results/BAM/"
	String VCF_path = "~/Singlepy/Results/VCF/"
	String input_path = "~/Singlepy/input"
	Int threads
}	

 call QualityCheck {
	input:
		sample_name = sample_name,
		r1fastq = r1fastq,
		r2fastq = r2fastq,	
	}
	
 call AlignBWA { 
	input:
		sample_name = sample_name,
		r1fastq = r1fastq,
		r2fastq = r2fastq,
		ref_fasta = ref_fasta,
		ref_fasta_amb = ref_fasta_amb,
		ref_fasta_sa = ref_fasta_sa,
		ref_fasta_bwt = ref_fasta_bwt,
		ref_fasta_ann = ref_fasta_ann,
		ref_fasta_pac = ref_fasta_pac,
	}
	
 call Markduplicates {
	input:
		sample_name = sample_name,
		outbam = AlignBWA.outbam	
	}
		
 call FixReadGroup {	
	input:
		sample_name = sample_name,
		out_dup_bam = Markduplicates.out_dup_bam
	}
		
call BuildBamIndex {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		bai_bam = FixReadGroup.out_fix	
	}
		
		
 call DepthOfCoverage {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_gene_list = ref_gene_list,
		interval_bed = interval_bed,
		input_bam = FixReadGroup.out_fix,
		depth_output = sample_name
	}
							
 call GATK_HaplotypeCaller {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		input_bam = FixReadGroup.out_fix,
		input_bam_index = BuildBamIndex.out_index,
                GATK_out = GATK_out,
		make_gvcf = make_gvcf,
                make_bamout = make_bamout,
		gatk_path = gatk_path
	}
		
 call select as selectSNP {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "SNP"
		F_vcf_SNP = GATK_HaplotypeCaller.GATK_out
	}
				
 call select as selectINDEL {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "INDEL"
		F_vcf_INDEL = GATK_HaplotypeCaller.GATK_out			
	}
	
 call VariantFiltration {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		input_vcf = GATK_HaplotypeCaller.GATK_out
	}
		
 call VariantRcalibrator_SNP {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		HapMap = HapMap
		omni = Omni
		1000G = 1000G
		dbSNP_hg19 = dbSNP_hg19
		VR_input_SNP = selectSNP.F_vcf_SNP
	}
		
 call VariantRecalibrator_INDEL {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		HapMap = HapMap
		omni = Omni
		1000G = 1000G
		dbSNP_hg19 = dbSNP_hg19
		VR_input_INDEL = selectINDEL.F_vcf_INDEL
	}
	
 call ApplyVQSR_SNP {
	input:
		sample_name = sample_name,
		VQSR_SNP =  selectSNP.F_vcf_SNP
		Recal_SNP = VariantRcalibrator_SNP.out_recal
	}
	
 Call ApplyVQSR_INDEL {
	input:
		sample_name = sample_name,
		VQRS_INDEL = selectINDEL.F_vcf_INDEL
		Recal_INDEL = VariantRecalibrator_INDEL.out_recal
	}
	
 call DeepCall {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		deep_input = FixReadGroup.out_dup_bam	
	}

 call select as selectSNP_deep {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "SNP"
		F_vcf_SNP_deep = DeepCall.sample_vcf
	}
				
 call select as selectINDEL_deep {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "INDEL"
		F_vcf_INDEL_deep = DeepCall.sample_vcf			
	}
	
call Bcftools_merge_SNP {
	input:
		sample_name = sample_name
		GATK_SNP = ApplyVQSR_SNP.output
		Deep_SNP = selectSNP_deep.output
	}

call Bcftools_merge_INDEl {
	input:
		sample_name = sample_name
		GATK_INDEL = ApplyVQSR_INDEL.output
		Deep_INDEL = selectINDEL_deep.output
	}

# Task 1 QC 
task QualityCheck {
    input {
	string sample_name
	File r1fastq
	File r2fastq
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
	Int threads = 6
}
	command <<<
	   set -x
	   # Check to see if input files are present
           [ -s ${r1fastq} ] || echo "Input Read1 Fastq is Empty" >> ${Failure_Logs}
           [ -s ${r2fastq} ] || echo "Input Read2 Fastq is Empty" >> ${Failure_Logs} 
	   
	   StartTime=`date +%s`
	   fastqc -t ${threads} ${r1fastq}  ${r2fastq} -o ${QC_path}
	   EndTime=`date +%s`
	   
	   echo -e "${sample_name} ran Quality check step for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
	   if [ $? -ne 0 ]; then
           echo "${sample_name} has failed at the FASTQ/BAM File Quality Control Step" >> ${Exit_Code}
        fi
           [ ! -d ${QC_path} ] && echo "QC directory has not been created" >> ${Failure_Logs}
      
	>>>

        runtime {
        continueOnReturnCode: true
   }	
}

# Task 2 this will do the alignment, Sort and convert sam to bam 
# Using BWA mem algorithm 

task AlignBWA {
   input {
	String sample_name
	File r1fastq
	File r2fastq
	File ref_fasta
	File ref_fasta_amb
	File ref_fasta_sa
	File ref_fasta_bwt
	File ref_fasta_ann
	File ref_fasta_pac
	Float mem_size_gb = 12
	Int threads
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		bwa mem \
		-t ${threads} \
		-R "@RG\tID:${sample_name}\tSM:${sample_name}" \
		${ref_path_hg19}/${ref_fasta} \
		${input_path}/${r1fastq} ${input_path}/${r2fastq} | \
		java -Xmx10g -jar ${picard_path}/picard.jar \
		SortSam \
		VALIDATION_STRINGENCY=SILENT \
		O= ${BAM_path}/~{sample_name}_hg19_sorted.bam \
		SORT_ORDER=coordinate
		EndTime=`date +%s`
		
		echo -e "${sample_name} Ran BWA Mem and Picard Sort for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		[[ ! -f ${BAM_path}/${sample_name}_hg19_sorted.bam ]] && echo -e "Aligned bam not created" >> ${Failure_Logs}
	>>>
	
	runtime {
		cpu: threads
		memory: "~{mem_size_gb} GiB"
	}
	
	output {
		
		File outbam = "${sample_name}_hg19_sorted.bam"
	}
	
}
	

# Task 3
## This task remove duplicates in BAM file

task Markduplicates {
	input {
	String sample_name
	File outbam
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		[ -s ${BAM_path}/${sample_name}_hg19_sorted.bam} ] || echo "Aligned sorted Bam File is Empty" >> ${Failure_Logs}
		
		StartTime=`date +%s`
		java -jar ${picard_path}/picard.jar \
		MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		AS=true \
		REMOVE_DUPLICATES=true \
		I= ${BAM_path}/~{outbam} \
		O= ${BAM_path}/~{sample_name}.hg19_sorted_mkdup.bam \
		M= ~{sample_name}.metrics
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran Picard Mark Duplicates for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${BAM_path}/~{sample_name}.hg19_sorted_mkdup.bam ]] && echo -e "Markduplicated bam not created" >> ${Failure_Logs}
	>>>
	
	output {
		
		File out_dup_bam = "~{sample_name}.hg19_sorted_mkdup.bam"
		File out_metrics = "~{sample_name}.metrics"
	}
	
}

# Task 4
## Add or fix Read Groups 

task FixReadGroup {
	input {
	String sample_name
	File out_dup_bam
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		[ -s ${hg19_sorted_mkdup.bam} ] || echo "Markduplicated Bam File is Empty" >> ${Failure_Logs}
		
		StartTime=`date +%s`
		java -jar ${picard_path}/picard.jar \
		AddOrReplaceReadGroups \
		VALIDATION_STRINGENCY=LENIENT \
		I= ${BAM_path}/~{out_dup_bam} \
		O= ${BAM_path}/~{sample_name}_hg19.bam \
		RGID=SRR"~{sample_name}" \
		RGLB=SRR"~{sample_name}" \
		RGPL=illumina \
		RGPU=SRR"~{sample_name}" \
		RGSM=SRR"~{sample_name}" \
		CREATE_INDEX=true
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran Picard AddOrReplaceReadGroups for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${BAM_path}/~{sample_name}_hg19.bam ]] && echo -e "ReadGroup Fixed bam not created" >> ${Failure_Logs}
	>>>
	
	output {
	
		File out_fix = "~{sample_name}_hg19.bam"
		
	}
	
}

# Task 5
## Build BAM file index 

task BuildBamIndex {
	input {
	String sample_name
	File out_fix
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		[ -s ${hg19.bam} ] || echo "ReadGroup fixed Bam File is Empty" >> ${Failure_Logs}
		StartTime=`date +%s`
		samtools index ${BAM_path}/~{out_fix}
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran Samtools index for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${BAM_path}/~{sample_name}_hg19.bam.bai ]] && echo -e "Bam Index File not created" >> ${Failure_Logs}
	
	>>>
	
	output {
	
		File out_index = "~{sample_name}_hg19.bam.bai"
		
	}
	
}

# Task 6
## Run DepthofCoverage for BAM file QC

task  DepthOfCoverage {
	input {
	String sample_name
	File out_fix
	File All_gene_hg19
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
	
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
		DepthOfCoverage \
		-R ${ref_path_hg19}/${ref_fasta} \
		-O ${BAM_path}/${sample_name} \
		-I ${BAM_path}/~{out_fix} \
		-gene-list ${ref_path_hg19}/${All_genes_hg19} \
		-L ${Input_files_target_bed}/truseq-dna-exome-targeted-regions-manifest-v1-2.bed

		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran DepthOfCoverage for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
	
	>>>
	
	output {
	
		File depth_output = "~{sample_name}"
		
	}
	
}

# Task 7
## Run GATK_HaplotypeCaller for variant calling

task GATK_HaplotypeCaller {
	input {
	String sample_name
	File out_fix
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
		HaplotypeCaller \
		-R ${ref_path_hg19}/${ref_fasta} \
		-I ${BAM_path}/~{out_fix} \
		-O ${VCF_path}/~{sample_name}_GATK.vcf.gz
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran GATK HaplotypeCaller for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/~{sample_name}_GATK.vcf.gz ]] && echo -e "GATK KaplotypeCaller VCF File not created" >> ${Failure_Logs}
	>>>
	
	output {
	
		File GATK_out = "~{sample_name}_GATK.vcf.gz"
		
	}
	
}

# Task 9
## Seperate SNPs from INDEL variants

task selectSNP {
	input {
	String sample_name
	File GATK_out
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		vcftools --gzvcf ${VCF_path}/~{GATK_out} \
		--remove-indels --recode --recode-INFO-all \
		--out ${VCF_path}/${sample_name}_GATK_SNP && mv ${VCF_path}/${sample_name}_GATK_SNP.recode.vcf ${VCF_path}/${sample_name}_GATK_SNP.vcf
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet SNPs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/~{sample_name}_GATK.vcf.gz ]] && echo -e "vcftools SNP sepeartion was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File GATK_SNP_out = "~{sample_name}_GATK_SNP.vcf"
		
	}
	
}

# Task 10
## Seperate INDELs from GATK VCF file

task selectINDEl {
	input {
	String sample_name
	File GATK_out
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}

	command <<<
		set -x
		
		StartTime=`date +%s`
		vcftools --gzvcf ${VCF_path}/~{GATK_out} \
		--remove-indels --recode --recode-INFO-all \
		--out ${VCF_path}/${sample_name}_GATK_INDEL && mv ${VCF_path}/${sample_name}_GATK_INDEL.recode.vcf ${VCF_path}/${sample_name}_GATK_INDEL.vcf
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/${sample_name}_GATK_INDEL.vcf ]] && echo -e "vcftools INDEL sepeartion was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File GATK_INDEL_out = "~{sample_name}_GATK_INDEL.vcf"
		
	}
	
}

# Task 11
## SNP Variants Hard Filtration Step

task VariantFiltration_SNP {
	input {
	String sample_name
	File GATK_SNP_out
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
		VariantFiltration \
		-R ${ref_path_hg19}/${ref_fasta} \
		-V ${VCF_path}/~{GATK_SNP_out} \
		-O ${VCF_path}/~{sample_name}_GATK_SNP_FP.vcf \
		--filter-expression "QD < 2.00" \
		--filter-name "QDlessthan2" --filter-expression "FS > 60.000" \
		--filter-name "FSgreaterthan60" --filter-expression "SOR > 3.000" \
		--filter-name "SORgreaterthan3" --filter-expression "MQRankSum<-2.50" \
		--filter-name "MQRankSumlessthanminus2.5" --filter-expression "MQRankSum > 2.50" \
		--filter-name "MQRankSumgreaaterthan2.5" --filter-expression "ReadPosRankSum < -1.00" \
		--filter-name "ReadPosRankSumlessthanminus1" --filter-expression "ReadPosRankSum > 3.50" \
		--filter-name "ReadPosRankSumgreaterthan3.5" --filter-expression "MQ < 60.00" \
		--filter-name "MQlessthan60"
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Varinat filterations using GATK for SNPs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/${sample_name}_GATK_SNP_FP.vcf ]] && echo -e "GATK Variant Hard Filteration for SNPs was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File GATK_SNP_FP = "~{sample_name}_GATK_SNP_FP.vcf"
		
	}
	
}

# Task 12
## INDEL Variants Hard Filtration Step

task VariantFiltration_INDEL {
	input {
	String sample_name
	File GATK_INDEL_out
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
		VariantFiltration \
		-R ${ref_path_hg19}/${ref_fasta} \
		-V ${VCF_path}/~{GATK_INDEL_out} \
		-O ${VCF_path}/~{sample_name}_GATK_INDEL_FP.vcf \
		--filter-expression "QD < 2.00" --filter-name "QDlessthan2" \
		--filter-expression "FS > 60.000" --filter-name "FSgreaterthan60" \
		--filter-expression "SOR > 3.000" --filter-name "SORgreaterthan3" \
		--filter-expression "ReadPosRankSum < -1.00" --filter-name "ReadPosRankSumlessthanminus1" \
		--filter-expression "ReadPosRankSum > 3.5" --filter-name "ReadPosRankSumgreaterthan3.5" \
		--filter-expression "MQ < 60.00" --filter-name "MQlessthan60"

		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Varinat filterations using GATK for INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/${sample_name}_GATK_INDEL_FP.vcf ]] && echo -e "GATK Variant Hard Filteration For INDELs was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File GATK_INDEL_FP = "~{sample_name}_GATK_INDEL_FP.vcf"
		
	}
	
}

# Task 13
## INDEL Variant Recalibration Step

task VariantRcalibrator_INDEL {
	input {
	String sample_name
	File mills_resource_vcf
        File dbsnp_vcf
        File mills_resource_vcf_index
        File dbsnp_vcf_index
	File GATK_INDEL_FP
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar VariantRecalibrator \
    		-V ${VCF_path}/~{GATK_INDEL_FP} \
    		--trust-all-polymorphic \
    		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
		-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
		-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    		-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
   		-mode INDEL \
    		--max-gaussians 4 \
   	 	-resource:mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
   	 	-resource:dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_vcf} \
   	 	-O ${VCF_path}/~{sample_name}_indels.recal \
   	 	--tranches-file ${VCF_path}/~{sample_name}_INDEL.tranches
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran VariantRecalibrator from GATK for INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/~{sample_name}_indels.recal ]] && echo -e "GATK Variant Recalibration For INDELs was not found" >> ${Failure_Logs}
	>>>

	output {
	
		File Recalib_INDEL = "~{sample_name}_GATK_INDEL_FP.vcf"
		File tranched_indel = "~{sample_name}_INDEL.tranches"
	}
	
}

# Task 14
## SNP Variant Recalibration step

task VariantRcalibrator_SNP {
	input {
	String sample_name
	File GATK_SNP_FP
	File hapmap_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf
        File omni_resource_vcf_index
        File thousand_G_resource_vcf
        File thousand_G_resource_vcf_index
        File dbsnp_vcf
        File dbsnp_vcf_index
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar VariantRecalibrator \
    		-V ${VCF_path}/~{GATK_SNP_FP} \
    		--trust-all-polymorphic \
    		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    		-mode SNP \
    		--max-gaussians 6 \
    		-resource:hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
    		-resource:omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
    		-resource:1000G,known=false,training=true,truth=false,prior=10:${thousand_G_resource_vcf} \
    		-resource:dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_vcf} \
    		-O ${VCF_path}/~{sample_name}_snp.recal \
    		--tranches-file ${VCF_path}/~{sample_name}_snps.tranches
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran VariantRecalibrator from GATK for SNPs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${VCF_path}/~{sample_name}_snp.recal ]] && echo -e "GATK Variant Recalibration For SNPs was not found" >> ${Failure_Logs}

	>>>
	
	output {
	
		File Recalib_SNP = "~{sample_name}_GATK_INDEL_FP.vcf"
		File tranches_snp = "~{sample_name}_snps.tranches"
		
	}
	
}
	
# Task 15
## ApplyVQSR step for SNP	

task ApplyVQSR_SNP {
	input {
	String sample_name
	File GATK_SNP_FP
	File Recalib_SNP
	File tranches_snp
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}	
	command <<<
		set -x
		
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
    		ApplyVQSR \
    		-V  ${VCF_path}/~{GATK_SNP_FP} \
    		--recal-file ${VCF_path}/${Recalib_SNP} \
    		--tranches-file ${VCF_path}/${tranches_snp} \
    		--truth-sensitivity-filter-level 99.7 \
    		--create-output-variant-index true \
    		-mode SNP \
    		-O ${VCF_path}/${sample_name}_SNP_recalibrated.vcf.gz 
		
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran ApplyVQSR from GATK for SNPs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${sample_name}_SNP_recalibrated.vcf.gz ]] && echo -e "GATK ApplyVQSR For SNPs was not found" >> ${Failure_Logs}
	>>>

	output {
	
		File VQSR_SNP = "~{sample_name}_SNP_recalibrated.vcf.gz"
		
	}
	
}

# Task 16
## ApplyVQSR step for INDEL

task ApplyVQSR_INDEL {
	input {
	String sample_name
	File GATK_INDEL_FP
	File Recalib_INDEL
	File tranched_indel
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		java -jar ${gatk_path}/gatk-package-4.2.4.0-local.jar \
   		ApplyVQSR \
    		-V ${VCF_path}/~{GATK_INDEL_FP} \
    		--recal-file ${VCF_path}/${Recalib_INDEL} \
    		--tranches-file ${VCF_path}/{tranched_indel} \
    		--truth-sensitivity-filter-level 99.7 \
    		--create-output-variant-index true \
    		-mode INDEL \
    		-O  ${VCF_path}/${sample_name}_indel_recalibrated.vcf.gz
		
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran ApplyVQSR from GATK for INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${sample_name}_indel_recalibrated.vcf.gz ]] && echo -e "GATK ApplyVQSR For INDELs was not found" >> ${Failure_Logs}
	>>>
	
	output {
	
		File VQSR_INDEL = "~{sample_name}_indel_recalibrated.vcf.gz"
		
	}
	
}

# Task 17
## Variant calling with DeepVariant 

task DeepCall {
	input {	
	String sample_name
	File out_fix
	File ref_fasta
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
	String mode_exome = "WES"
	String mode_genome = "WGS"
	Int threads
}
	command <<<
		set -x
		StartTime=`date +%s`
		singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  		docker://google/deepvariant:1.3.1 \
  		/opt/deepvariant/bin/run_deepvariant \
  		--model_type= ${mode_exome} \ 
  		--ref=  ${ref_path_hg19}/${ref_fasta} \
  		--reads= ${BAM_path}/~{out_fix} \
  		--output_vcf= ${VCF_path}/~{sample_name}_deepVariant.vcf.gz \
  		--output_gvcf= ${VCF_path}/~{sample_name}_deepVariant.g.vcf.gz \
  		--intermediate_results_dir ${VCF_path} \ 
  		--num_shards= ${threads} \ 
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran DeepVariant for variant calling for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${sample_name}_deepVariant.vcf.gz ]] && echo -e "Variant calling with DeepVariant was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File Deep_vcf = "~{sample_name}_GATK_INDEL_FP.vcf"
		
	}
	
}

# Task 18
## Seperating SNPS from DeepCall output

task selectSNP_deep {
	input {
	String sample_name
	File Deep_vcf
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		vcftools --gzvcf ${VCF_path}/~{Deep_vcf} \
		--remove-indels --recode --recode-INFO-all \
		--out ${VCF_path}/${sample_name}_Deep_SNP && mv ${sample_name}_Deep_SNP.recode.vcf ${sample_name}_Deep_SNP.vcf
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet SNPs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${sample_name}_Deep_SNP.vcf ]] && echo -e "DeepVariant SNP variant seperation was not successful" >> ${Failure_Logs}
	>>>
	
	output {
	
		File deep_SNP = "~{sample_name}_Deep_SNP.vcf"
		File deep_SNP_record = "~{sample_name}_Deep_SNP.recode.vcf"
		
	}
	
}
				
# Task 19
## Seperating INDELs from DeepCall output		
task selectINDEL_deep {
	input {
	String sample_name
	File Deep_vcf
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		
		StartTime=`date +%s`
		vcftools --gzvcf ${VCF_path}/~{Deep_vcf} \
		--remove-indels --recode --recode-INFO-all \
		--out ${VCF_path}/${sample_name}_Deep_Indel && mv ${sample_name}_Deep_Indel.recode.vcf ${sample_name}_Deep_Indel.vcf
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		[[ ! -f ${sample_name}_Deep_Indel.vcf ]] && echo -e "DeepVariant Indel variant seperation was not successful" >> ${Failure_Logs}
		
	>>>
	
	output {
	
		File deep_Indel = "~{sample_name}_Deep_Indel.vcf"
		File deep_Indel_record = "~{sample_name}_Deep_Indel.recode.vcf"
	}
	
}
		
# Task 20
## Merging GATK and DeepVariant SNP outputs

task Bcftools_merge_SNP {
	input {
	String sample_name
	File deep_SNP
	File VQSR_SNP
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		bcftools merge ${VCF_path}/${VQSR_SNP} ${VCF_path}/${deep_SNP} > ${VCF_path}/${sample_name}_merged_snp.vcf
		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
	>>>
	
	output {
	
		File GATK_INDEL_out = "~{sample_name}_GATK_INDEL.vcf"
		
	}
	
}
		
		
task Bcftools_merge_INDEl {
	input {
	String sample_name
	File deep_Indel
	File VQSR_INDEL
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		StartTime=`date +%s`
		bcftools merge file1.vcf.gz fle2.vcf.gz file3.vcf.gz > out.vcf
				EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran vcftools to separet INDELs for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
	>>>
	
	output {
	
		File GATK_INDEL_out = "~{sample_name}_GATK_INDEL.vcf"
		
	}
	
}

