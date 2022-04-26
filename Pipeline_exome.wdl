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
	File ref_gene_list
	File interval_bed
	File ref_dict
	File ref_index
	File dbSNO_vcf
	File dbSNP_vcf_index
	File All_gene
	File HapMap
	File Omni
	File 1000G
	File dbSNP_hg19
	Array[File] known_indels_vcf
	Array[File] known_indels_indices
	Boolean make_gvcf = true
	String gatk_path = ""
	String QC_path = ""
	String BAM_path = ""
	String Var_path = ""
	String Tools_path = ""
	String Input_files_all_gene = ""
	String Input_files_target_bed = ""
	
	
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
	
 call SortSam {
	input : 
		sample_name = sample_name,
		insam = AlignBWA.outsam	
	}
		
 call Markduplicates {
	input:
		sample_name = sample_name,
		outbam = SortSam.outbam	
	}
		
 call FixReadGroup {	
	input:
		sample_name = sample_name,
		out_dup_bam = Markduplicates.outbam
	}
		
call BuildBamIndex {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		bai_bam = FixReadGroup.out_dup_bam	
	}
		
		
 call DepthOfCoverage {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_gene_list = ref_gene_list,
		interval_bed = interval_bed,
		input_bam = FixReadGroup.out_dup_bam,
		depth_output = sample_name
	}
		
 call Qualimap {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_gene_list = ref_gene_list,
		interval_bed = interval_bed,
		input_bam = FixReadGroup.out_dup_bam,	
	}
							
 call GATK_HaplotypeCaller {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		input_bam = FixReadGroup.out_dup_bam
		input_bam_index = BuildBamIndex.bai_bam
                output_filename = output_filename,
		make_gvcf = make_gvcf,
                make_bamout = make_bamout,
		gatk_path = gatk_path
	}
		
 call VariantFiltration {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		input_vcf = GATK_HaplotypeCaller.sample_vcf
	}
	
 call select as selectSNP {
	input:
		sample_name = sample_name,
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "SNP"
		F_vcf_SNP = VariantFiltration.sample_vcf
	}
				
 call select as selectINDEL {
	input:
		sample_name = sample_name
		ref_fasta = ref_fasta,
		ref_dict = ref_dict,
		ref_index = ref_index,
		type = "INDEL"
		F_vcf_INDEL = VariantFiltration.sample_vcf			
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

# Tasks #
## Align reads to the reference genome hg19/hg38
## using BWA algorithm

task QualityCheck {
    input {
	string sample_name
	File r1fastq
	File r2fastq
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command {
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

      [ ! -d ${QC_path} ] && echo "FASTQC directory has not been created" >> ${Failure_Logs}
      
	}

        runtime {
      continueOnReturnCode: true
   }	
}

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
		bwa mem -t ${threads} -R "@RG\tID:${sample_name}\tSM:${sample_name}" ${ref_fasta} ${r1fastq} ${r2fastq} -o ${BAM_path}/${sample_name}.hg19-bwamem.sam
		EndTime=`date +%s`
		
		echo -e "${sample_name} Ran BWA Mem for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		[[ ! -f ${sample_name}.aligned.bam ]] && echo -e "Aligned bam not created" >> ${Failure_Logs}
	>>>
	
	runtime {
		cpu: threads
		memory: "~{mem_size_gb} GiB"
	}
	
	output {
		
		File outsam = "${sample_name}.hg19-bwamem.sam"
	}
	
}
	
# Task2
## This task will sort and convert sam to bam file and index the BAM

task SortSam {
	input {
	String sample_name
	File insam
	String Exit_Code
	String Failure_Logs
	String dollar = "$"
}
	command <<<
		set -x
		# To see whether the sam files is created or not!
		[ -s ${hg19-bwamem.sam} ] || echo "Aligned Sam File is Empty" >> ${Failure_Logs}
		StartTime=`date +%s`
		
		java -Xmx10g -jar ${Tools_path}/picard.jar \
		SortSam \
		VALIDATION_STRINGENCY=SILENT \
		I= ${BAM_path}/~{insam} \
		O= ${BAM_path}/~{sample_name}.hg19-bwamem.sorted.bam \
		SORT_ORDER=coordinate
		
		EndTime=`date +%s`
		echo "${sample_name} Ran Samtools SortSam for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
		
		 [ ! -f ${sample_name}.hg19-bwamem.sorted.bam ] && echo "Aligned sorted bam not created" >> ${Failure_Logs}
		
	>>>
	
	output {
	
		File outbam = "~{sample_name}.hg19-bwamem.sorted.bam"
	}
	
}

# Task3
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
		[ -s ${hg19-bwamem.sorted.bam} ] || echo "Aligned sorted Bam File is Empty" >> ${Failure_Logs}
		
		StartTime=`date +%s`
		java -jar ${Tools_path}/picard.jar \
		MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		AS=true \
		REMOVE_DUPLICATES=true \
		I= ${BAM_path}/~{outbam} \
		O= ${BAM_path}/~{sample_name}.hg19-bwamem.sorted.mkdup.bam \
		M= ~{sample_name}.metrics
		EndTime=`date +%s`
		# How long dose it take
		echo "${sample_name} Ran Picard Mark Duplicates for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
	>>>
	
	output {
		
		File out_dup_bam = "~{sample_name}.hg19-bwamem.sorted.mkdup.bam"
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
		[ -s ${hg19-bwamem.sorted.mkdup.bam} ] || echo "Markduplicated Bam File is Empty" >> ${Failure_Logs}
		StartTime=`date +%s`
		java -jar ${Tools_path}/picard.jar \
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
	>>>
	
	output {
	
		File out_fix = "~{sample_name}_hg19.bam"
		
	}
	
}

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
	
	>>>
	
	output {
	
		File out_index = "~{sample_name}_hg19.bai"
		
	}
	
}
		
task  DepthOfCoverage {
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
		java -jar ${Tools_path}/gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar \
		DepthOfCoverage \
		-R ${ref_fasta} \
		-O ${BAM_path}/${sample_name} \
		-I ${BAM_path}/~{out_fix} \
		-gene-list ${Input_files_all_gene}/All_gene.refseq \
		-L ${Input_files_target_bed}/truseq-dna-exome-targeted-regions-manifest-v1-2.bed

		EndTime=`date +%s`
		
		# How long dose it take
		echo "${sample_name} Ran DepthOfCoverage for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
	
	>>>
	
	output {
	
		File out_index = "~{sample_name}"
		
	}
	
}






