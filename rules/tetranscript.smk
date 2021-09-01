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
    params:
        format          =       config["TEtranscript"]["format"],
        stranded        =       config["TEtranscript"]["stranded"],
        project         =       config["TEtranscript"]["project"]
    shell:
        "TEtranscripts -t {input.treatment} -c {input.control} --GTF {input.genic_gtf} --TE {input.TE_gtf} --format {params.format} --stranded {params.stranded} --project {params.project}"

