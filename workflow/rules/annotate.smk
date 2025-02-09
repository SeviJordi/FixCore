rule annotate:
    threads: config["PROKKA"]["N_CORES"]
    conda: "../envs/prokka.yaml" # TODO: create env
    params:
        outdir = PROKKA_DIR
    input:
        genome = Path(config["GENOMES_DIR"])/"{genome_name}.fasta"
    output:
        gff = PROKKA_DIR/"{genome_name}/{genome_name}.gff"
    shell:
        """
        prokka --cpus {threads} \
            --outdir {params.outdir}/{wildcards.genome_name} \
            --prefix {wildcards.genome_name} {input.genome}
        """
