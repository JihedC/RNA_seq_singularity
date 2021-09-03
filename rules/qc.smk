################## Rules used for QC ##################


rule fastp:
    input:
        get_fastq
    output:
        fq1  = temp(WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz"),
        fq2  = temp(WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"),
        html = WORKING_DIR + "fastp/{sample}_fastp.html",
        json = WORKING_DIR + "fastp/{sample}_fastp.json"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        in_and_out_files =  get_trim_names,
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    resources: cpus=10
    shell:
        "touch {output.fq2};\
        fastp --thread {threads}  --html {output.html} --json {output.json} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        {params.in_and_out_files} \
        2>{log}"


rule multiqc:
    input:
        expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "

rule generate_mapping_summary:
    input:
        expand(RESULT_DIR + "star/{sample}_Log.final.out", sample = SAMPLES)
    output:
        RESULT_DIR + "mapping_summary.csv"
    message:
        "Concatenating STAR mapping report and generating .csv mapping summary."
    params:
        directory_with_mapping_reports = RESULT_DIR + "star/" 
    shell:
        "python scripts/generate_mapping_summary.py {params.directory_with_mapping_reports} {output}"
