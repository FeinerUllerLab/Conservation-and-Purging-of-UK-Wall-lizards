#!/bin/bash

# Objective: Prepare the reference genome for the downstream analyses

module load samtools
samtools faidx genome.fna

module load bwa
bwa index genome.fna

PICARD_HOME=/pdc/software/eb/software/picard/2.25.5
java -Xmx16g -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=genome.fna O=genome.dict
