#parameters for the star alignment are based on : https://static-content.springer.com/esm/art%3A10.1186%2Fs13100-019-0192-1/MediaObjects/13100_2019_192_MOESM5_ESM.pdf

rule mapping:
	input:
		fastq	=	get_fastq,
		#indexdone="{indexDirectory}/index.DONE".format(indexDirectory=config["reference"]["index"]),
		#annotation= "{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	params:
		STAR="--alignEndsType EndToEnd --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonical",
		prefix="results/mapped/{samples}/{samples}",
		genome="/exports/humgen/jihed/TEtranscript/mm10"
	threads:
		32
	log:
		"results/log/star/{samples}.log"
	conda:
		"../envs/star.yaml"
	output:
		"results/mapped/{samples}/{samples}Aligned.sortedByCoord.out.bam"
	shell:
		"""
		STAR --readFilesIn {input.fastq} --runMode alignReads {params.STAR} --runThreadN {threads} --genomeDir {params.genome} --outFileNamePrefix {params.prefix}
		"""

#STAR="--outSAMtype BAM Unsorted --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000",

#"STAR --runMode alignReads {params.STAR} --outFileNamePrefix {params.prefix} --runThreadN {threads} --sjdbGTFfile {input.annotation} --genomeDir {params.genome} --readFilesIn {input.fastq} 2>{log}"

rule sort_bam:
	input:"results/mapped/{samples}/{samples}Aligned.sortedByCoord.out.bam",
	output: "results/mapped/{samples}/{samples}.sorted.bam"
	log:
		"results/log/sort/{samples}.log"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools sort -n -o {output} {input}"

rule download_gtf_gene:
	output:
		"annotation/annotation_gene.gtf"
	params:
		gtfFile=config["GENOME_ZIP_GTF_URL"]
	log:
		"results/log/annotation/download_gene_annotation.log"
	shell:
		"curl {params.gtfFile} | gunzip -c > {output}"

rule download_gtf_repeat:
	output:
		"annotation/annotation_repeat.gtf"
	params:
		gtfFile=config["REPEAT_GTF_URL"]
	log:
		"results/log/annotation/download_repeat_annotation.log"
	shell:
		"curl {params.gtfFile} | gunzip -c > {output}"

rule deduplicate:
		input:
			bam="results/mapped/{samples}/{samples}Aligned.sortedByCoord.out.bam"
		output:
			dedup="results/mapped/{samples}/{samples}.dedup.bam",
			stats="results/mapped/{samples}/{samples}.dedup.stats"
		conda:
			"../envs/picard.yaml"
		shell:
			"picard MarkDuplicates I={input.bam} O={output.dedup} M={output.stats} REMOVE_DUPLICATES=true"
#picard MarkDuplicates I=results/mapped/WT1/WT1Aligned.out.bam O=results/mapped/WT1/WT1Aligned.out.bam METRICS_FILE=test.picard.txt REMOVE_DUPLICATES=true

# java -Xmx10g -jar /u/project/jacobsen/resources/scripts_and_pipelines/scripts/MarkDuplicates.jar
# I="$read1" O="$outdir/STAR/${name}_dedup.bam"
# METRICS_FILE="$outdir/STAR/other_files_and_logs/${name}_MarkDuplicates_metrics.txt"
# REMOVE_DUPLICATES=true > "$outdir/STAR/other_files_and_logs/${name}_MarkDuplicates_log.txt" 2>&1
rule htseq_count:
		input:
			lambda wildcards: expand("results/mapped/{samples}/{samples}.dedup.bam", samples = SAMPLES)
		output:
			"results/counts/htseq_count.txt"
		params:
			gtf_repeat	=	"GRCm38_GENCODE_rmsk_TE.gtf"
		log:
			"results/log/htseqcount/results/htseqcount_log.txt"
		shell:
			"htseq-count --format=bam --idattr=transcript_id{input} {params.gtf_repeat} > {output}"

# rule TEtranscripts:
# 		input:
# 			sorted		=  	"results/mapped/{samples}/{samples}.sorted.bam"
# 			gtf_gene	=	"gencode.vM20.annotation.gtf"
# 			gtf_repeat	=	"GRCm38_GENCODE_rmsk_TE.gtf"
# 		output:
# 			"results/te_local/{samples}.cntTable"
# 		log:
# 			"results/log/te_local/{samples}.log"
# 		shell:
# 			"TElocal --sortByPos -b {input.sorted} --GTF {input.gtf_gene} --TE {input.gtf_repeat} --project {output}"
