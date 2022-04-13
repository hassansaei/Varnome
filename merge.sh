#!/bin/bash

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

	do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_R1.fastq.gz

	echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_R2.fastq.gz

done;

echo "Alignment with BWA-MEM"
for fname in *_R1.fastq.gz; do base=${fname%_R1*}; bwa mem -t 8 ~/quickstart-ref/BWA_index/hg19_ucsc.fa  "${base}_R1.fastq.gz"  "${base}_R2.fastq.gz" >"${base}.sam"; done

echo "SAM to BAM"
for f  in *.sam; do  samtools view -@ 10 -bS $f  > ${f/.sam/.bam} ; done

echo "Sort SAM"
for f  in *.bam; do java -Xmx10g -jar ~/bin/picard.jar SortSam VALIDATION_STRINGENCY=SILENT I=$f O=${f/.bam/_Sort.bam} SORT_ORDER=coordinate; done

echo "MarkDuplicate"
for f in *_Sort.bam; do java -jar ~/bin/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT AS=true REMOVE_DUPLICATES=true I=$f O=${f/_Sort.bam/_mkdup.bam} M=${f/_Sort.bam/.metrics}; done



