rule concat_alignments:
    threads: 8
    conda:
        "../envs/phylo.yaml"
    params:
        temp_dir = TEMP_DIR
    input:
        fastas = expand(PATHCURATED/"{gene_name}.mafft.evalmsa.trimmal.fasta", gene_name=iter_gene_names())
    output:
        concatenated = OUTDIR/f"{PREFIX}.concatenated.fixcore.fasta",
        partitions = OUTDIR/f"{PREFIX}.concatenated.fixcore.partitions"
    log:
        LOGDIR/"concat_alignments.log"
    shell:
        """      
        exec >{log}
        exec 2>&1

        # Move non empty files to temp dir

        for file in {input.fastas}; do
            nlines=$(wc -l $file | cut -d" " -f1)
            nseqs=$(grep -c ">" $file)

            if [[ $nlines -gt $nseqs ]]; then
                mv $file {params.temp_dir}
            fi
            
        done

        AMAS.py concat -f fasta -d dna -i "{params.temp_dir}"/* -c 8 -t {output.concatenated} -p {output.partitions}

        rm -r {params.temp_dir}
        """

rule create_vcf:
    conda:
        "../envs/phylo.yaml"
    threads: 1
    input:
        alignment = OUTDIR/f"{PREFIX}.concatenated.fixcore.fasta"
    output:
        vcf = PATHPHYLO/f"{PREFIX}.vcf",
    log:
        LOGDIR/"create_vcf.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        snp-sites -v {input.alignment} > {output.vcf} 
        """


rule get_SNPs_alignment:
    conda:
        "../envs/phylo.yaml"
    threads: 1
    input:
        alignment = OUTDIR/f"{PREFIX}.concatenated.fixcore.fasta"
    output:
        snps = PATHPHYLO/f"{PREFIX}.SNPs.fasta",
    log:
        LOGDIR/"get_SNPs_alignment.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        snp-sites {input.alignment} > {output.snps}
        """

rule iqtree:
    threads: 8
    conda:
        "../envs/phylo.yaml"
    params:
        prefix = f"{PATHPHYLO}/{PREFIX}",
        model = config["EVO_MODEL"],
        bootstrap_replicates = config["BOOTSTRAP"]
    input:
        snps = PATHPHYLO/f"{PREFIX}.SNPs.fasta",
        alignment = OUTDIR/f"{PREFIX}.concatenated.fixcore.fasta"
    output:
        treefile = PATHPHYLO/f"{PREFIX}.treefile",
        constantvar = PATHPHYLO/f"{PREFIX}.fconst.txt"
    log:
        LOGDIR/"iqtree.log"
    shell:
        """
        exec >{log}
        exec 2>&1
        
        fconstsvar=$(snp-sites -C  {input.alignment})
        echo $fconstsvar > {output.constantvar}

        iqtree2 -s {input.snps} \
            -m {params.model} \
            -B {params.bootstrap_replicates} \
            -nt {threads} \
            --redo -fconst $fconstsvar \
            --prefix {params.prefix}
        """

        