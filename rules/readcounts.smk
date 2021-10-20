################## Rules used for read counting ##################


rule htseq_count:
    input:
        lambda wildcards: expand(RESULT_DIR + "star/{sample}_Aligned.out.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "htseq/htseq_count_TE.txt"
    params:
        TE_gtf	=	WORKING_DIR + "TE_repeat_masker.gtf",
        header  =   "Gene\\\t"+"\\\t".join(samples)
    log:
        RESULT_DIR + "log/htseq_count/htseq_count.log"
    message:
        "Producing TE count table with HTSEQ-count"    
    shell:
        """
        echo {params.header} > {output}
        htseq-count --format=bam --idattr=transcript_id {input} {params.TE_gtf} >> {output}
        """

rule featurecount_TE:
    input:
        bams = expand(RESULT_DIR + "star/{sample}_Aligned.out.bam", sample = SAMPLES),
        TE_gtf	=	WORKING_DIR + "TE_repeat_masker.gtf"
    output:
        RESULT_DIR + "featureCounts/featureCounts_TE.txt"
    log:
        RESULT_DIR + "log/featureCounts/featureCounts_TE.log"
    message:
        "Producing TE count table with featureCounts"
    threads: 10    
    shell:
        "featureCounts -T {threads} -t exon -g transcript_id -F 'gtf' -a {input.TE_gtf} -o {output} {input.bams}"


rule featurecount_genes:
    input:
        bams = expand(RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        gene_gtf	=	WORKING_DIR + "annotation.gtf"
    output:
        RESULT_DIR + "featureCounts/genes/featureCounts_genes.txt"
    log:
        RESULT_DIR + "log/featureCounts/featureCounts_genes.log"
    message:
        "Producing Genes count table with featureCounts"
    threads: 10    
    shell:
        "featureCounts -T {threads} -F 'gtf' -a {input.gene_gtf} -o {output} {input.bams}"       

rule createCountsPerRepetitiveRegions:
	input:
		bamFiles    =   expand(RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
		annotation  =   annotation + "mm10.rm.bed.gz"
	output:
		RESULT_DIR + "Global_TE.countsPerRepetitiveRegions.csv"
	params:
		header="chr\\\tstart\\\tend\\\tID\\\t\\\tsize\\\tstrand\\\t"+"\\\t".join(SAMPLES)
	shell:
		"""
		echo {params.header}>{output}
		bedtools multicov -bams {input.bamFiles} -bed {input.annotation}>> {output}
		"""
