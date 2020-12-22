rule mapping:
	input:
		fastq	=	getFastq,
		indexdone="{indexDirectory}/index.DONE".format(indexDirectory=config["reference"]["index"]),
		annotation= "{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	params:
		STAR="--outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType bedGraph --outWigNorm RPM --readFilesCommand zcat  --quantMode GeneCounts --twopassMode Basic",
		prefix="analysis/mappedDataDir/{sample}/{sample}",
		genome=config["reference"]["index"]
	threads:
		config["threads"]
	log:
		"results/log/star/{samples}.log"
	conda:
		"../envs/star.yaml"
	output:

	shell:
		"""
		STAR --readFilesIn {input} \
			--runThreadN 4 --outSAMtype BAM Unsorted --runMode alignReads \
		  	--outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 \
  			--outMultimapperOrder Random --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd \
  			--alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 \
  			--alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 \
  			--seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 --genomeDir mm10/ --outFileNamePrefix WT
		"""


rule sort_bam:
	input:
	output:
	log:
	shell:
	""

rule te_transcript:
		input:
		output:
		log:
		conda:
		shell:
		""
