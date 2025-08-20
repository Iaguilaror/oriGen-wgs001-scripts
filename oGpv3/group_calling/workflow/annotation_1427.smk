import pandas as pd
configfile: "../config/config.yaml"

dataset = pd.read_csv(config["dataset"], sep="\t")
RUNIDS = dataset["runid"]
CODEIDS = dataset["codeid"]
CRAM_VERSIONS = dataset["cram_version"]
CHROMS = config["chroms"]
cram_version_lookup = dataset.set_index(["runid", "codeid"])["cram_version"]

ANNOVAR_DBS = {
    "refGene": "g",
    "cytoBand": "r",
    "dbnsfp31a_interpro": "f",
    "gnomad40_genome": "f",
    "clinvar_20221231": "f",
    "dbnsfp42a": "f",
    "avsnp150": "f",
    "intervar_20180118": "f",
    "revel": "f",
    "mcap": "f"
}

localrules: all, generate_sams_file

wildcard_constraints:
    varnum = "(_var[0-9]+)?"

def additional_threads(wildcards, threads):
    return max(threads - 1, 1)

def threads_minus_two(wildcards, threads):
    return max(threads - 2, 1)

rule all:
    input:
        expand("../results/RELEASE_PHASE1/FASE1_oriGen-1427_annoSNV_{chrom}{ext}",
            chrom=CHROMS, ext=[".vcf.gz", ".vcf.gz.tbi"])

rule annotate:
    input:
        ref = config["ref_genome"],
        dbdir = "../resources/annovar_dbs/",
        vcf = "../results/RELEASE_PHASE1/FASE1_oriGen-1427_SNV_{chrom}.vcf.gz"
    output:
        avinput = temp("../work/snv-anno_{chrom}.avinput"),
        vcf = temp("../work/snv-anno_{chrom}.hg38_multianno.vcf"),
        txt = temp("../work/snv-anno_{chrom}.hg38_multianno.txt"),
    benchmark: "../benchmarks/annotate/annotate-{chrom}.tsv"
    log: "../logs/annotate/annotate-{chrom}.log"
    threads: config["anno_threads"]
    params:
        dbs = ",".join(ANNOVAR_DBS.keys()),
        ops = ",".join(ANNOVAR_DBS.values()),
        athreads = additional_threads
    shell: '''
            outfile={output.avinput}
            bcftools view -G -Ov {input.vcf} | \\
            table_annovar.pl \\
                --outfile "${{outfile%.*}}" \\
                --buildver hg38 \\
                --protocol {params.dbs} \\
                --operation {params.ops} \\
                --vcfinput \\
                --remove \\
                --thread {params.athreads} \\
                --maxgenethread {params.athreads} \\
                --polish \\
                /dev/stdin {input.dbdir} > {log} 2>&1
    '''

rule reheader_anno_snv:
    input:
        bcf = "../results/RELEASE_PHASE1/FASE1_oriGen-1427_SNV_{chrom}.vcf.gz",
        bcfi = "../results/RELEASE_PHASE1/FASE1_oriGen-1427_SNV_{chrom}.vcf.gz.tbi",
        colmap = "../resources/column_remapping.txt",
        header = "../resources/vcf_info_header.txt",
        vcf = "../work/snv-anno_{chrom}.hg38_multianno.vcf"
    output:
        vcf = "../results/RELEASE_PHASE1/FASE1_oriGen-1427_annoSNV_{chrom}.vcf.gz",
        vcfi = "../results/RELEASE_PHASE1/FASE1_oriGen-1427_annoSNV_{chrom}.vcf.gz.tbi"
    benchmark: "../benchmarks/reheader_anno/reheader_anno-{chrom}.tsv"
    log: "../logs/reheader_anno/reheader_anno-{chrom}.log"
    threads: config["anno_threads"]
    priority: 50
    params:
        athreads = additional_threads
    shell: '''
            (
                bcftools view -Gh {input.bcf} | \\
                perl -ane '
                    if($f < 2) {{
                        if(/^##FILTER=/) {{ $f = 1 }}
                        elsif($f == 1) {{
                            $f = 2;
                            open($h, "<", "{input.header}") or die;
                            while($l = readline($h)) {{print "$l"}};
                            close($h);
                        }}
                    }}
                    print';
                
                cat {input.vcf} | \\
                perl -ane '
                    BEGIN {{
                        open($h, "<", "{input.colmap}") or die;
                        while($l = readline($h)) {{
                            chomp $l;
                            @a = split("\\t", $l);
                            $m{{$a[0]}} = $a[1];
                        }}
                        close($h);
                    }}
                    for $k (keys(%m)) {{
                        s/(?<=;)\\Q${{k}}\\E(?=\\=)/${{m{{$k}}}}/g
                    }}
                    print';
            ) | \\
                bcftools view \\
                    --write-index=tbi \\
                    --threads {params.athreads} \\
                    -Oz \\
                    -o {output.vcf} >> {log} 2>&1
    '''