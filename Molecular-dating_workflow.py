###########################################################
# 2023-03-29
# First version of molecular dating workflow
###########################################################

from Bio import SeqIO


REFERENCE = "NC_007596.2"
FASTA = "DVP-mito-NDexclded-3old"
PRIORS = "Priors_DPVetal-mitogenomes_5oldest_5youngest"
SAMPLES = []

with open(FASTA+'.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            taxa = str(line).strip('>').strip('\n')
            taxa_name = taxa.split('_')[0]
            taxa_date = taxa.split('_')[-1]
            if taxa_date == "ND":
                SAMPLES.append(taxa_name)


rule all:
    input:
        expand(FASTA + "_{sample}.log", sample=SAMPLES),
        expand(FASTA + "_{sample}.trees", sample=SAMPLES),
        expand(FASTA + "_{sample}.out", sample=SAMPLES),
        expand(FASTA + "_{sample}.xml", sample=SAMPLES),



rule creat_target_seq_for_liftoff:
    input:
        fasta = FASTA + ".fasta",
        ref = "Liftoff/ref_" + REFERENCE + ".fasta"
    output:
        target = "Liftoff/target_" + REFERENCE + ".fasta"
    run:
        fa_dict = SeqIO.to_dict(SeqIO.parse(FASTA + ".fasta", 'fasta'))
        target_dict = {}
        for k, s in fa_dict.items():
            if REFERENCE in k:
                target_dict[k] = s
            else:
                short = REFERENCE.split(".")[0]
                if short in k:
                    target_dict[REFERENCE] = s
        with open("Liftoff/target_" + REFERENCE + ".fasta", 'w') as handle:
            SeqIO.write(target_dict.values(), handle, 'fasta')



rule adjust_annotation:
    input:
        fasta = FASTA + ".fasta",
        ref = "Liftoff/ref_" + REFERENCE + ".fasta",
        target = "Liftoff/target_" + REFERENCE + ".fasta",
        gff = "Liftoff/" + REFERENCE + ".gff3",
        partitions = "Liftoff/partitions"
    output:
        "Liftoff/" + REFERENCE + ".liftoff.gff3"
    log:
        "Liftoff/" + REFERENCE + ".liftoff.log"
    shell:
        "liftoff {input.target} {input.ref} -g {input.gff} -o {output} -f {input.partitions} > {log} 2>&1"


rule generate_xml:
    input:
        fasta = FASTA + ".fasta",
        template = "data/template.xml",
        gff = "Liftoff/" + REFERENCE + ".liftoff.gff3",
        priors = "data/" + PRIORS + ".csv"
    output:
        xml = FASTA + "_{sample}.xml"
    params:
        fa = FASTA + ".fasta"
    log:
        "generatorLog/" + FASTA + "_{sample}.XMLgenerator.log"
    shell:
        "python XMLgenerator_Ver1.1.py \
        -f {params.fa} \
        -t {input.template} \
        -g {input.gff} \
        -p {input.priors} \
        -m 10000000 -l 1000 -i True \
        > {log} 2>&1"
 

rule run_beast:
    input:
        xml = FASTA + "_{sample}.xml"
    output:
        multiext(FASTA + "_{sample}", ".log", ".trees")
    log:
        "runLog/" + FASTA + "_{sample}.beast_run.log"
    shell:
        "beast -beagle_cpu -beagle_sse -beagle_double \
        -threads {threads} {input.xml} > {log} 2>&1"



rule tree_annotator:
    input:
        FASTA + "_{sample}.trees"
    output:
        FASTA + "_{sample}.out"
    log:
        FASTA + "_{sample}.treeAnnotator.log"
    shell:
        "treeannotator -burnin 1000000 {input} {output}"

