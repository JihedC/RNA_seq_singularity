################## Rules used for TE transcript ##################

rule TEtranscripts:
    input:
        treatment       =       expand(RESULT_DIR + "star/{treatment}_Aligned.out.bam", treatment = CASES),
        control         =       expand(RESULT_DIR + "star/{control}_Aligned.out.bam", control = CONTROLS),
        genic_gtf       =       WORKING_DIR + "annotation.gtf",
        TE_gtf          =       WORKING_DIR + "TE_repeat_masker.gtf"
    output:
        RESULT_DIR + "TEtranscript/TEtranscript_out.cntTable"
    log:
        RESULT_DIR + "log/TEtranscript/tetranscript.log"
	singularity:'docker:/mhammelllab/tetranscripts:latest'        
    params:
        format          =       config["TEtranscript"]["format"],
        stranded        =       config["TEtranscript"]["stranded"],
        project         =       config["TEtranscript"]["project"]
    shell:
        "TEtranscripts -t {input.treatment} -c {input.control} --GTF {input.genic_gtf} --TE {input.TE_gtf} --format {params.format} --stranded {params.stranded} --project {params.project}"


rule TElocal:
    input:
        bam		    =  	RESULT_DIR + "star/{sample}_Aligned.out.bam",
        genic_gtf	=	WORKING_DIR + "annotation.gtf",
        TE_gtf	    =	WORKING_DIR + "TE_prebuilt_index.locInd"
    output:
        RESULT_DIR + "TElocal/{sample}.cntTable.cntTable"
    params:
        project     =   RESULT_DIR + "TElocal/{sample}.cntTable",
        stranded    =   config["TElocal"]["stranded"],
    shell:
        "TElocal -b {input.bam} --GTF {input.genic_gtf} --TE {input.TE_gtf} --project {params.project} --stranded {params.stranded}"
