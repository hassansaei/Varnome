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
	File bed_library_prep_kit
	File All_gene
	Array[File] known_indels_vcf
	Array[File] known_indels_indices
}	
	call qualityCheck {
		input:
			sample_name = sample_name,
			r1fastq = r1fastq,
			r2fastq = r2fastq,
		}
	call alignBWA { 
	
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
	
	call sortSam {
		input : 
			sample_name = sample_name,
			insam = alignBWA.outsam
			
		}
		
	call Markduplicates {
		input:
			sample_name = sample_name,
			outbam = sortSam.outbam
			
		}
		
	call DepthofCoverage {
		input:
			sample_name = sample_name,
			ref_fasta = ref_fasta,
			All_gene = All_gene,
			bed_library_prep_kit = bed_library_prep_kit,
			indepth = Markduplicates.outbam
		
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
			FixedBAM = FixReadGroup.out_dup_bam
			
		}
		
	call HaplotypeCaller {
		
		input:
			sample_name = sample_name,
			ref_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			BamFile = 
			BamFile_index =
			
		}
		
	call select as selectSNP {
		
		input:
			sample_name = sample_name,
			ref_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			type = "SNP"
			rawVCF=HaplotypeCaller.sample_vcf
		}
		
	call select as selectINDEL {
	
		input:
			sample_name = sample_name
			ref_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			type = "INDEL"
			rawVCF=HaplotypeCaller.sample_vcf
			
		}
	
	call HardFilterSNP {
	
		input:
			sample_name = sample_name,
			ref_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			rawSNP= selectSNP.rawvcf
			
		}
		
	call HardFilterINDEL {
	
		input:
		
			sample_name = sample_name,
			ref_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			rawSNP= selectINDEL.rawvcf
			
		}
		
	call DeepVariant {
	
		input:
			sample_name = sample_name
			ef_fasta = ref_fasta,
			ref_dict = ref_dict,
			ref_index = ref_index,
			
}



# Tasks #
## Align reads to the reference genome hg19/hg38
## using BWA algorithm

task align {
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
	Int threads
}
	
	command { 
		bwa mem -t ${threads} -R "@RG\tID:A8100\tSM:A8100" ${ref_fasta} ${r1fastq} ${r2fastq} -o ${sample_name}.hg19-bwamem.sam 
	}
	
	runtime {
		cpu: threads
		requested_memory_mb: 16000
	}
	
	output {
		
		File outsam = "${sample_name}.hg19-bwamem.sam"
	}
	
}
	
	
# Task2
## This task will sort and convert sam to bam file and index the BAM

task sortSam {

	input {
	
	String sample_name
	File insam
}
	command <<<
		java -Xmx10g -jar ~/bin/picard.jar \
		SortSam \
		VALIDATION_STRINGENCY=SILENT \
		I= ~{insam} \
		O= ~{sample_name}.hg19-bwamem.sorted.bam \
		SORT_ORDER=coordinate
		
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
}
	command <<<
		java -jar ~/bin/picard.jar \
		MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		AS=true \
		REMOVE_DUPLICATES=true \
		I=~{outbam} \
		O=~{sample_name}.hg19-bwamem.sorted.mkdup.bam \
		M=~{sample_name}.metrics
		
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
}
	command <<<
		java -jar ~/bin/picard.jar \
		AddOrReplaceReadGroups \
		VALIDATION_STRINGENCY=LENIENT \
		I=~{out_dup_bam} \
		O=~{sample_name}.hg19-bwamem.bam \
		RGID=SRR"~{sample_name}" \
		RGLB=SRR"~{sample_name}" \
		RGPL=illumina \
		RGPU=SRR"~{sample_name}" \
		RGSM=SRR"~{sample_name}" \
		CREATE_INDEX=true
		
	>>>
	
	output {
	
		File out_fix = "~{sample_name}.hg19-bwamem.bam"
		
	}
	
}


