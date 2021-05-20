from glob import glob 
configfile: "config.yml"
FASTQ_DIR= config["FASTQ_DIR"]
SAMPLES = [re.search(FASTQ_DIR+r"/(.+)_1.fastq.gz",i).group(1) for i in glob(os.path.join(FASTQ_DIR,"*_1.fastq.gz"))]
OUT_DIR = config["OUT_DIR"]

rule pangolin:
    conda: 
        "../envs/illumina/environment.yml",
    
    input:
        expand("{out_dir}/consensus/{sample}.primertrimmed.consensus.fa", sample=SAMPLES, out_dir=[OUT_DIR])
    
    output:
        fasta=OUT_DIR +"/annotation/all_consensus_seqs.fa",
        csv=OUT_DIR +"/annotation/pangolin_lineage_report.csv"

    shell:
        """
        cat {input} > {output.fasta}
        pangolin {output.fasta} --outfile {output.csv}
        """

rule anno_CorGat:
    conda: 
        "../envs/illumina/environment.yml",

    # label 'error_ignore'
    input:
        consensus=OUT_DIR+"/consensus/{sample}.primertrimmed.consensus.fa"

    output:
        align=OUT_DIR+"/annotation/{sample}_align.csv",
        anno=OUT_DIR+"/annotation/{sample}_anno.csv"
    
    params:
        corgatPath=config['corgatPath']+'/',
        corgatFna=config['corgatFna'],
        corgatConf=config['corgatConf'],
        test_conf=OUT_DIR+"/annotation/test.conf",
        align=os.path.join(config['corgatPath'], 'align.pl'),
        annotate=os.path.join(config['corgatPath'], 'annotate.pl'),
        tmp=OUT_DIR+"/tmp/"

     
    shell:
        """
        if [ ! -f {params.test_conf} ];then
            awk -v path={params.corgatPath} '{{full_path=path$2; print $1,full_path}}' {params.corgatConf} > {params.test_conf}
        fi
        mkdir -p {params.tmp}{wildcards.sample}
        cd {params.tmp}{wildcards.sample}
        {params.align} --multi  {input.consensus} --refile {params.corgatFna}  --out {output.align}
        {params.annotate} --in {output.align} --conf {params.test_conf} --out {output.anno}
        """