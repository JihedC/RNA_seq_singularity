################## Rules used to download genome and GTF annotations ##################


rule download_genome:
	params: 
		fasta   =   config["GENOME_ZIP_FASTA_URL"]
	output:
		WORKING_DIR + "reference.fa"
	log:
		"results/log/download_genomefile.log"    
	shell:
		"curl {params.fasta} | gunzip -c > {output}"


rule download_gtf_gene:
	output:
		WORKING_DIR + "annotation.gtf"
	params:
		gtfFile=config["GENOME_ZIP_GTF_URL"]
	log:
		"results/log/download_gene_annotation.log"
	shell:
		"curl {params.gtfFile} | gunzip -c > {output}"

rule download_TE_gene:
	output:
		WORKING_DIR + "TE_repeat_masker.gtf"
	params:
		gtfFile=config["REPEAT_GTF_URL"]
	log:
		"results/log/download_TE_repeat_masker.log"
	shell:
		"curl {params.gtfFile} | gunzip -c > {output}"
