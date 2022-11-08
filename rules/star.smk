################## Rules used to map the data with STAR ##################

rule star_index:
    input:
        fasta = WORKING_DIR + "reference.fa",
        gtf =   WORKING_DIR + "annotation.gtf"
    output:
         genome_index = [WORKING_DIR + "genome/" + f for f in ["chrLength.txt","chrNameLength.txt","chrName.txt","chrStart.txt","Genome","genomeParameters.txt","SA","SAindex"]]
    message:
        "generating STAR genome index"
    params:
        genome_dir = WORKING_DIR + "genome/",
        sjdb_overhang = config["star_index"]["sjdbOverhang"],
    threads: 10
    resources: mem_mb=100000
    singularity:'docker:/mgibio/star:latest'        
    shell:
        "mkdir -p {params.genome_dir}; " # if directory not created STAR will ask for it
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.genome_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "

rule map_to_genome_using_STAR:
    input:
        genome_index = rules.star_index.output,
        forward_read = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        reverse_read = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        gtf =   WORKING_DIR + "annotation.gtf"
    output:
        RESULT_DIR + "star/{sample}_Aligned.out.bam",
        RESULT_DIR + "star/{sample}_Log.final.out"
    message:
        "mapping {wildcards.sample} reads to genome"
    log:
        RESULT_DIR + "log/star/{sample}.log"
    benchmark:
        RESULT_DIR + "benchmark/star_{sample}_unsorted.benchmark.txt"
    singularity:'docker:/mgibio/star:latest'        
    params:
        sample_name           =  "{sample}",
        star_input_file_names =  get_star_names,
        prefix                =  RESULT_DIR + "star/{sample}_",
        maxmismatches         =  config["star"]["mismatches"],
        unmapped              =  config["star"]["unmapped"]   ,
        multimappers          =  config["star"]["multimappers"],
        matchNminoverLread    =  config["star"]["matchminoverlengthread"],
        outSamType            =  config["star"]["samtype"],
        outSAMattributes      =  config["star"]["samattributes"],
        intronmax             =  config["star"]["intronmax"],
        matesgap              =  config["star"]["matesgap"],
        genome_index          =  WORKING_DIR + "genome/",
        sjdbOverhang          =  config["star"]["sjdbOverhang"],
        winAnchorMultimapNmax =  config["star"]["winAnchorMultimapNmax"]    
    threads: 10
    resources: cpus=10
    shell:
        "STAR --runThreadN 12 --genomeDir {params.genome_index} --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.sjdbOverhang} --readFilesIn {params.star_input_file_names} \
        --readFilesCommand zcat --winAnchorMultimapNmax {params.winAnchorMultimapNmax} \
        --outFilterMultimapNmax {params.multimappers} \
        --outReadsUnmapped {params.unmapped} \
        --outFileNamePrefix {params.prefix} --outSAMtype {params.outSamType} "

rule map_to_genome_using_STAR_sorted:
    input:
        genome_index = rules.star_index.output,
        forward_read = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        reverse_read = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        gtf =   WORKING_DIR + "annotation.gtf"
    output:
        RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam",
        RESULT_DIR + "sorted_star/{sample}_Log.final.out"
    message:
        "mapping {wildcards.sample} reads to genome using David's parameters"
    log:
        RESULT_DIR + "log/star/{sample}.log"
    benchmark:
        RESULT_DIR + "benchmark/star_{sample}_sorted.benchmark.txt"        
    params:
        sample_name           =  "{sample}",
        star_input_file_names =  get_star_names,
        prefix                =  RESULT_DIR + "sorted_star/{sample}_",
        maxmismatches         =  config["star"]["mismatches"],
        unmapped              =  config["star"]["unmapped"]   ,
        multimappers          =  config["star"]["multimappers"],
        matchNminoverLread    =  config["star"]["matchminoverlengthread"],
        outSamType_sorted     =  config["star"]["samtype_sorted"],
        outSAMattributes      =  config["star"]["samattributes"],
        intronmax             =  config["star"]["intronmax"],
        matesgap              =  config["star"]["matesgap"],
        genome_index          =  WORKING_DIR + "genome/",
        sjdbOverhang          =  config["star"]["sjdbOverhang"],
        winAnchorMultimapNmax =  config["star"]["winAnchorMultimapNmax"]    
    threads: 10
    resources: cpus=10
    singularity:'docker:/mgibio/star:latest'
    shell:
        "STAR --runThreadN 12 \
        --genomeDir {params.genome_index} \
        --sjdbGTFfile {input.gtf} \
        --outFileNamePrefix {params.prefix} \
        --readFilesIn {params.star_input_file_names} \
        --readFilesCommand zcat \
        --outSAMtype {params.outSamType_sorted} \
        --bamRemoveDuplicatesType UniqueIdentical \
        --outWigType bedGraph --outWigNorm RPM \
        --quantMode GeneCounts --twopassMode Basic"
        


rule index_bam:
    input:
        RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    log:
        RESULT_DIR + "log/sort/{sample}.log"
    singularity:'docker:/biocontainers/samtools:v1.9-4-deb_cv1'
    shell:
        "samtools index {input} 2>{log}"

rule bamcoverage:
    input:
        bam     =   RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam",
        bai     =   RESULT_DIR + "sorted_star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        bigwig  =   RESULT_DIR + "bigwig/{sample}_rpkm.bw"
    message:
        "Create genome coverage tracks"
    benchmark:
        RESULT_DIR + "benchmark/bamcoverage_{sample}.benchmark.txt"
    params:
         binsize                =   config["bamcoverage"]["binsize"],
         normalizeUsing         =   config["bamcoverage"]["normalizeUsing"],
         effectiveGenomeSize    =   config["bamcoverage"]["effectiveGenomeSize"],
         smoothLength           =   config["bamcoverage"]["smoothLength"]       
    log:
        RESULT_DIR + "log/bamcoverage/{sample}.log"
    singularity:'docker://stjudecloud/deeptools:branch-chipseq-1.0.2'     
    shell:
        "bamCoverage -b {input.bam} --binSize {params.binsize} --effectiveGenomeSize {params.effectiveGenomeSize} --normalizeUsing {params.normalizeUsing} --smoothLength {params.smoothLength} -o {output.bigwig} 2>{log}"