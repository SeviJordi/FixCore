rule concat_alignments:
    threads: 8
    conda:
        "../envs/phylo.yaml"
    input:
        fastas = expand(PATHCURATED/"{gene_name}.mafft.evalmsa.trimmal.fasta", gene_name=get_gene_names())
    output:
        concatenated = OUTDIR/"{PREFIX}.concatenated.fixcore.fasta"
        partitions = OUTDIR/"{PREFIX}.concatenated.fixcore.partitions"
    shell:
        """
        AMAS.py concat -f fasta -d dna -i {input.fastas} -c 8 -t {output.concatenated} -p {output.partitions}
        """

rule create_vcf:
    conda:
        "../envs/phylo.yaml"
    input:
        alignment = OUTDIR/"{PREFIX}.concatenated.fixcore.fasta"
    output:
        vcf = PATHPHYLO/"{PREFIX}.vcf"
    shell:
        """
        snp-sites -v {input.alignment} > {output.vcf}
        """


rule get_SNPs_alignment:
    conda:
        "../envs/phylo.yaml"
    input:
        alignment = OUTDIR/"{PREFIX}.concatenated.fixcore.fasta"
    output:
        snps = PATHPHYLO/"{PREFIX}.SNPs.fasta"
    shell:
        """
        snp-sites {input.alignment} > {output.snps}
        """

rule iqtree:
    threads: 8
    conda:
        "../envs/phylo.yaml"
    params:
        prefix = f"{PATHPHYLO}/{PREFIX}"
        model = config["EVO_MODEL"]
        bootstrap_replicates = config["BOOTSTRAP"]
    input:
        snps = ATHPHYLO/"{PREFIX}.SNPs.fasta"
    output:
        treefile = PATHPHYLO/"{PREFIX}.treefile"
        constantvar = PATHPHYLO/"{PREFIX}.fconst.txt"
    shell:
        """
        fconstsvar=$(snp-sites -C  $PATHTREE/"$PREFIX"_fixcore.fasta)
        echo $fconstsvar > {output.constantvar}

        iqtree2 -s $PATHTREE/"$PREFIX"_fixcore.SNPs.fasta \
            -m {params.model} \
            -B {params.bootstrap_replicates} \
            -nt {threads} \
            --redo -fconst $fconstsvar \
            --prefix {params.prefix}
        """

        