#parameters for the star alignment are based on : https://static-content.springer.com/esm/art%3A10.1186%2Fs13100-019-0192-1/MediaObjects/13100_2019_192_MOESM5_ESM.pdf

rule mapping:
	input:
		fastq	=	get_fastq,
		#indexdone="{indexDirectory}/index.DONE".format(indexDirectory=config["reference"]["index"]),
		#annotation= "{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	params:
		STAR="--outSAMtype BAM Unsorted --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000",
		prefix="results/mapped/{samples}/{samples}",
		genome="/exports/humgen/jihed/TEtranscript/mm10"
	threads:
		8
	log:
		"results/log/star/{samples}.log"
	conda:
		"../envs/star.yaml"
	output:
		"results/mapped/{samples}/{samples}Aligned.out.bam"
	shell:
		"""
		STAR --readFilesIn {input.fastq} --runMode alignReads {params.STAR} --runThreadN {threads} --genomeDir {params.genome} --outFileNamePrefix {params.prefix}
		"""
#"STAR --runMode alignReads {params.STAR} --outFileNamePrefix {params.prefix} --runThreadN {threads} --sjdbGTFfile {input.annotation} --genomeDir {params.genome} --readFilesIn {input.fastq} 2>{log}"

rule sort_bam:
	input:"results/mapped/{samples}/{samples}.Aligned.out.bam"
	output: "results/mapped/{samples}/{samples}.sorted.bam"
	log:
		"results/log/sort/{samples}.log"
	conda:
		"../conda/samtools.yaml"
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

rule te_local:
		input:
			sorted		=  	"results/mapped/{samples}/{samples}.sorted.bam"
			gtf_gene	=	"gencode.vM20.annotation.gtf"
			gtf_repeat	=	"mm10_rmsk_TE.gtf.locInd"
		output:
			"results/te_local/{sample}.cntTable"
		log:
			"results/log/te_local/{sample}.log"
		shell:
			"TElocal --sortByPos -b {input.sorted} --GTF {input.gtf_gene} --TE {input.gtf_repeat} --project {output}"
