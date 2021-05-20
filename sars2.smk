configfile: "config.yml"
include: "rules/anno.smk"

from glob import glob 
import re
import os


# include: "rules/illumina.smk"
# include "rules/qc.smk"

def get_final_output():
    final_output = []
    final_output.extend(expand(OUT_DIR+"/annotation/{sample}_anno.csv", sample=SAMPLES)),
    

    return final_output

# Global variable 
GENOME = config["GENOME"]
GENOME_BWA_INDEX = GENOME + ".bwt"
GENOME_FAI_INDEX = GENOME + ".fai"
GENOME_NAME="NC_045512.2"
OUT_DIR = config["OUT_DIR"]
PRIMER_DIR = config["PRIMER_DIR"]
PRIMER_BED = os.path.join(PRIMER_DIR, "nCoV-2019.primer.bed")
# Get samples names from fastq 
FASTQ_DIR= config["FASTQ_DIR"]

SAMPLES = [re.search(FASTQ_DIR+r"/(.+)_1.fastq.gz",i).group(1) for i in glob(os.path.join(FASTQ_DIR,"*_1.fastq.gz"))]

print(SAMPLES)


rule root:
    input:
        get_final_output(),
        OUT_DIR +"/annotation/pangolin_lineage_report.csv",
        OUT_DIR + "/qc/quast/report.txt",
        OUT_DIR +"/qc/qc_plots/total_qc.csv",
        OUT_DIR + "/typing_summary.csv",
        OUT_DIR + "/variant_summary.csv",
        OUT_DIR +'/multiqc_report.html'


rule fastqc:
    conda: 
        "envs/illumina/environment.yml",
    input:
        R1=FASTQ_DIR + "/{sample}_1.fastq.gz",
        R2=FASTQ_DIR + "/{sample}_2.fastq.gz",
    output:
        html1=OUT_DIR + "/qc/fastqc/{sample}_1_fastqc.html",
        html2=OUT_DIR +"/qc/fastqc/{sample}_2_fastqc.html",
    params:
        out_dir=OUT_DIR + "/qc/fastqc/",
    shell:
        """
        fastqc -o {params.out_dir} -f fastq -q {input.R1} {input.R2}
        """  

rule readTrimming:
    conda: 
        "envs/illumina/environment.yml",

    input:
        R1=FASTQ_DIR + "/{sample}_1.fastq.gz",
        R2=FASTQ_DIR + "/{sample}_2.fastq.gz",

    output:
        OUT_DIR +"/trimmed/{sample}_1_val_1.fq.gz",
        OUT_DIR +"/trimmed/{sample}_2_val_2.fq.gz",
    log:
        OUT_DIR +"/logs/trimmed/{sample}_1.log",
        OUT_DIR +"/logs/trimmed/{sample}_2.log",
    params:
        out_dir=OUT_DIR +"/trimmed"
    shell:
        "trim_galore --paired {input.R1} {input.R2} -o {params.out_dir}"


rule index_genome:
    conda: "envs/illumina/environment.yml"
    input:
        GENOME
    output:
        GENOME_FAI_INDEX,
        GENOME_BWA_INDEX
    shell:
        "bwa index {input}; samtools faidx {input}"

rule bwa_align:
    conda: "envs/illumina/environment.yml"
    input:
        R1=OUT_DIR +"/trimmed/{sample}_1_val_1.fq.gz",
        R2=OUT_DIR +"/trimmed/{sample}_2_val_2.fq.gz",
        index=GENOME_BWA_INDEX

    output:
        OUT_DIR +"/mapped/{sample}.sam"
    log:
        OUT_DIR +"/logs/mapped/{sample}.log"
    params:
        genome = GENOME
    shell:
        "bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' {params.genome} {input.R1} {input.R2} > {output} 2> {log}"


rule sam2bam:
    conda: "envs/illumina/environment.yml"
    input:
        OUT_DIR +"/mapped/{sample}.sam"
    output:
        OUT_DIR +"/bams/{sample}.bam"
    shell:
        "samtools sort -O BAM {input} > {output}"


rule samIndex:
    conda: "envs/illumina/environment.yml"
    input:
        OUT_DIR +"/bams/{sample}.bam"
    output:
        OUT_DIR +"/bams/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule trimPrimerSequences:
    conda: 
        "envs/illumina/environment.yml",
    input:
        bam=OUT_DIR +"/bams/{sample}.bam",
        bai=OUT_DIR +"/bams/{sample}.bam.bai",
        bedfile=PRIMER_BED

    output:
        mapped=OUT_DIR +"/mapped/{sample}.mapped.bam",
        ptrim=OUT_DIR +"/mapped/{sample}.primertrimmed.sorted.bam",

    log:
        OUT_DIR +"/mapped/{sample}_trimPrimerSequences.log",

    params:
        ivarCmd = "ivar trim -e" if config["allowNoprimer"] else "ivar trim",
        illuminaKeepLen = config["illuminaKeepLen"],
        illuminaQualThreshold = config["illuminaQualThreshold"],
        cleanBamHeader = config["cleanBamHeader"]
    
    script:
        "bin/trimPrimerSequences.py"


rule callVariants:
    conda: 
        "envs/illumina/environment.yml",

    input:
        bam=OUT_DIR +"/mapped/{sample}.primertrimmed.sorted.bam",
        genome=GENOME

    output:
        OUT_DIR +"/variants/{sample}.variants.tsv",

    params:
        ivarMinDepth=config["ivarMinDepth"],
        ivarMinVariantQuality=config["ivarMinVariantQuality"],
        ivarMinFreqThreshold=config["ivarMinFreqThreshold"]

    shell:
        """
        samtools mpileup -A -d 0 --reference {input.genome} -B -Q 0 {input.bam} |\
        ivar variants -r {input.genome} -m {params.ivarMinDepth} -p {wildcards.sample}.variants -q {params.ivarMinVariantQuality} -t {params.ivarMinFreqThreshold} 
        mv  {wildcards.sample}.variants.tsv {output}
        """

rule makeConsensus:
    conda: 
        "envs/illumina/environment.yml",

    input:
        bam=OUT_DIR +"/mapped/{sample}.primertrimmed.sorted.bam",

    output:
        fasta=OUT_DIR +"/consensus/{sample}.primertrimmed.consensus.fa",
        quality=OUT_DIR +"/consensus/{sample}.primertrimmed.consensus.qual.txt"
        # quallity=directory("consensus")
    params:
        mpileupDepth=config["mpileupDepth"],
        ivarFreqThreshold=config["ivarFreqThreshold"],
        ivarMinDepth=config["ivarMinDepth"],


    shell:
        """
        samtools mpileup -aa -A -B -d {params.mpileupDepth} -Q0 {input.bam} | \
        ivar consensus -t {params.ivarFreqThreshold} -m {params.ivarMinDepth} \
        -n N -p {wildcards.sample}.primertrimmed.consensus.fa
        mv {wildcards.sample}.primertrimmed.consensus.fa {output.fasta}
        mv {wildcards.sample}.primertrimmed.consensus.qual.txt {output.quality}
        """

rule ivar_QUAST:
    conda: 
        "envs/illumina/environment.yml",

    input:
        consensus=expand("{out_dir}/consensus/{sample}.primertrimmed.consensus.fa", sample=SAMPLES, out_dir=[OUT_DIR]),
        # consensus=OUT_DIR +"/annotation/all_consensus_seqs.fa",
        genome=GENOME,
        

    output:
        out_dir=directory(OUT_DIR + "/qc/quast"),
        txt=OUT_DIR + "/qc/quast/report.txt"

    params:
        features = f"--features {config['gff']}" if config['gff'] else ""

    shell:
        """
        quast.py \\
            --output-dir {output.out_dir} \\
            -r {input.genome} \\
            {params.features} \\
            {input.consensus}
        """


rule makeQCCSV:
    conda: 
        "envs/illumina/environment.yml",

    input:
        bam=OUT_DIR +"/mapped/{sample}.primertrimmed.sorted.bam",
        consensus=OUT_DIR+"/consensus/{sample}.primertrimmed.consensus.fa",
        genome=GENOME


    output:
        csv=OUT_DIR +"/qc/qc_plots/{sample}.qc.csv",
        png=OUT_DIR +"/qc/qc_plots/{sample}.depth.png",
        

    params:
        qcSetting = "--illumina" if config["illumina"] else "--nanopore"

    shell:
        """
        ./bin/qc.py {params.qcSetting} --outfile {output.csv} --sample {wildcards.sample} --ref {input.genome} --bam {input.bam} --fasta {input.consensus}
        mv {wildcards.sample}.depth.png {output.png}
        """

rule writeQCSummaryCSV:
    conda: 
        "envs/illumina/environment.yml",
    input:
        csv=expand("{out_dir}/qc/qc_plots/{sample}.qc.csv", sample=SAMPLES, out_dir=[OUT_DIR])
    output:
        total_qc=OUT_DIR +"/qc/qc_plots/total_qc.csv"
    shell:
        "awk '(NR == 1) || (FNR > 1)' {input.csv} > {output.total_qc}"

# NR for total number of rows for all files in {input.csv}
# FNR for number of rows for each file in {input.csv}

rule typeVariants:
    conda: 
        "envs/illumina/environment.yml",

    input:
        genome=GENOME,
        variants=OUT_DIR +"/variants/{sample}.variants.tsv",
        gff=config['gff'],
        yaml=config['typing_yaml'],

    output:
        variants=OUT_DIR + "/typing/{sample}.variants.csv",
        typing=OUT_DIR + "/typing/{sample}.typing.csv",
        csq=OUT_DIR + "/typing/{sample}.csq.vcf"

    params:
        csqDpThreshold=config['csqDpThreshold'],
        csqAfThreshold=config['csqAfThreshold'],
        seq_type='-t' if config['illumina'] else '-v'

    shell:
        """
        ./bin/type_vcf.py \
        -i {wildcards.sample} \
        -y {input.yaml} \
        -ov {output.csq} \
        -ot {output.typing} \
        -os {output.variants} \
        -dp {params.csqDpThreshold} \
        -af {params.csqAfThreshold} \
        {params.seq_type} {input.variants} \
        {input.gff} {input.genome}
        """


rule mergeTypingCSVs:
    conda: 
        "envs/illumina/environment.yml",

    input:
        variants=expand("{out_dir}/typing/{sample}.variants.csv",sample=SAMPLES, out_dir=[OUT_DIR]),
        typing=expand("{out_dir}/typing/{sample}.typing.csv",sample=SAMPLES, out_dir=[OUT_DIR]),

    output:
        typing=OUT_DIR + "/typing_summary.csv",
        variants=OUT_DIR + "/variant_summary.csv"

    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input.variants} > {output.variants}
        awk '(NR == 1) || (FNR > 1)' {input.typing} > {output.typing}
        """


rule multiQC:
    conda: 
        "envs/illumina/environment.yml",
    
    input:
        fastqc=expand("{out_dir}/qc/fastqc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=[1,2], out_dir=[OUT_DIR]),
        mapped=expand("{out_dir}/mapped/{sample}_trimPrimerSequences.log", sample=SAMPLES, out_dir=OUT_DIR),
        quast=OUT_DIR + "/qc/quast/report.txt"
    
    output:
         OUT_DIR +'/multiqc_report.html'

    log:
        OUT_DIR + "/logs/multiqc.log"
     
    script:
        "bin/multiqc_wrapper.py"
