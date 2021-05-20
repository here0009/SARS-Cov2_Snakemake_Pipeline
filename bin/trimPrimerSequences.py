#! usr/bin/env
from snakemake.shell import shell

if snakemake.params.cleanBamHeader:
    shell(
        """
        samtools reheader --no-PG  -c 'sed "s/{snakemake.wildcards.sample}/sample/g"' {snakemake.input.bam} | \
        samtools view -F4 -o {snakemake.output.mapped}
   
        samtools index {snakemake.output.mapped}

        {snakemake.params.ivarCmd} -i {snakemake.output.mapped} -b {snakemake.input.bedfile} -m {snakemake.params.illuminaKeepLen} -q {snakemake.params.illuminaQualThreshold} -p ivar_{snakemake.wildcards.sample} > {snakemake.log}

        samtools reheader --no-PG  -c 'sed "s/{snakemake.wildcards.sample}/sample/g"' ivar_{snakemake.wildcards.sample}.bam | \
        samtools sort -o {snakemake.output.ptrim}
        rm -f ivar_{snakemake.wildcards.sample}.bam
        """
    )
else:
    shell(
        """
        samtools view -F4 -o {snakemake.output.mapped} {snakemake.input.bam}
        samtools index {snakemake.output.mapped}
        {snakemake.params.ivarCmd} -i {snakemake.output.mapped} -b {snakemake.input.bedfile} -m {snakemake.params.illuminaKeepLen} -q {snakemake.params.illuminaQualThreshold} -p ivar_{snakemake.wildcards.sample} > {snakemake.log}
        samtools sort -o {snakemake.output.ptrim} ivar_{snakemake.wildcards.sample}.bam
        rm -f ivar_{snakemake.wildcards.sample}.bam
        """
    )