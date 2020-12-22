################## Import libraries ##################

import pandas as pd
import os
import sys
from subprocess import call
import itertools
from snakemake.utils import R


################## Configuration file and PATHS##################

configfile: "config.yaml"

units = pd.read_table(config["units"], dtype=str).set_index(["bed"], drop=False)

BED = units.index.get_level_values('bed').unique().tolist()

units = pd.read_table(config["sample"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

################## DESIRED OUTPUT ##################

MAPPED  = ""
SORTED  = ""

################## RULE ALL ##################

rule all:
    input:
        MAPPED,
        SORTED

    message : "Analysis is complete!"
    shell:""


################## INCLUDE RULES ##################


include: "rules/te_transcript.smk"
