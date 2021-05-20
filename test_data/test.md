# Test the workflow

```bash
# convert nextflow pipeline to snakemake pipeline
# rename .nf to .smk
for f in *.nf; do 
    mv -- "$f" "${f%.nf}.smk"
done 

# run the pipeline
snakemake -s sars2.smk --use-conda -j 4
snakemake -s sars2.smk -n

trim_galore --paired SRR14476228_1.fastq.gz  SRR14476228_2.fastq.gz -o trimmed  --basename SRR14476228 > trim.log

SRR=SRR14476228
samtools reheader --no-PG  -c 'sed "s/SRR14476229/sample/g"' bams/SRR14476229.bam | samtools view -F4 -o mapped/SRR14476229.mapped.bam

samtools index mapped/SRR14476229.mapped.bam

ivar trim -e -i mapped/SRR14476229.mapped.bam -b /home/dlf/Code/SnakeMake/artic_ncov2019/ref/nCoV-2019/V3/nCoV-2019.primer.bed -m 30 -q 20 -p ivar_SRR14476229 > logs/primer_trime/SRR14476229_trimPrimerSequences.log

samtools reheader --no-PG  -c 'sed "s/SRR14476229/sample/g"' ivar_SRR14476229.bam |samtools sort -o mapped/SRR14476229.primertrimmed.sorted.bam

samtools sort -o mapped/SRR14476229.primertrimmed.sorted.bam' ivar_SRR14476229.bam'

awk -v path=/home/dlf/Code/SnakeMake/artic_ncov2019/CorGAT '{full_path=path\$2; print \$1,full_path}' /home/dlf/Code/SnakeMake/artic_ncov2019/CorGAT/corgat.conf > /home/dlf/Code/SnakeMake/artic_ncov2019/results/annotation/test.conf

 awk -v path=/home/dlf/Code/SnakeMake/artic_ncov2019/CorGAT/ '{full_path=path$2; print $1,full_path}' /home/dlf/Code/SnakeMake/artic_ncov2019/CorGAT/corgat.conf > /home/dlf/Code/SnakeMake/artic_ncov2019/results/annotation/test.conf


if [ ! -f /home/dlf/Code/SnakeMake/artic_ncov2019/results/annotation/test.conf ];then
    echo "file not exist"
fi


rule mergeTypingCSVs:

    conda: 
        "envs/illumina/environment.yml",

    input:
        variants=expand("{out_dir}/typing/{sample}.variants.csv",sample=SAMPLES, out_dir=[OUT_DIR]),
        typing=expand("{out_dir}/typing/{sample}.typing.csv",sample=SAMPLES, out_dir=[OUT_DIR]),

    output:
        typing=OUT_DIR + "typing_summary.csv",
        variants=OUT_DIR + "variant_summary.csv"

    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input.variants} > {output.variants}
        awk '(NR == 1) || (FNR > 1)' {input.typing} > {output.typing}
        """
