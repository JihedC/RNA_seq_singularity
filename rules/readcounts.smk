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

rule featurecount:
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