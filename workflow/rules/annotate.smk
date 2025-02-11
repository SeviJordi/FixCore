
rule annotate:
    threads: config["PROKKA"]["N_CORES"]
    conda: "../envs/prokka.yaml" # TODO: create env
    params:
        outdir = PROKKA_DIR
    input:
        genome = Path(config["GENOMES_DIR"])/"{genome_name}.fasta"
    output:
        gff = PROKKA_DIR/"{genome_name}/{genome_name}.gff"
    log:
        LOGDIR/"annotate"/"{genome_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        sed "s/ .*//g" {input.genome} > {input.genome}.renamed

        prokka --cpus {threads} \
            --outdir {params.outdir}/{wildcards.genome_name} \
            --prefix {wildcards.genome_name} {input.genome}.renamed --force
        
        rm {input.genome}.renamed
        """
