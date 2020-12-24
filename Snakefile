################## Import libraries ##################

import pandas as pd
import os
import sys
from subprocess import call
import itertools
from snakemake.utils import R


################## Configuration file and PATHS##################

configfile: "config.yaml"

#units = pd.read_table(config["units"], dtype=str).set_index(["bed"], drop=False)

#BED = units.index.get_level_values('bed').unique().tolist()

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.samples), ["fq1", "fq2"]].dropna()


##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"
################## DESIRED OUTPUT ##################

SORTED      = expand("results/mapped/{samples}/{samples}.sorted.bam", samples=SAMPLES)
TE_local    = expand("results/te_local/{samples}.cntTable", samples=SAMPLES)
################## RULE ALL ##################

rule all:
    input:
        SORTED

    message : "Analysis is complete!"
    shell:""


################## INCLUDE RULES ##################


include: "rules/te_transcript.smk"
