rule run_trimal:
    shadow:
        "shallow"
    conda:
        "../envs/trimal.yaml"
    threads: 1
    params:
        threshold = config["TRIMAL_PCT"]
    input:
        aligned = PATHCURATED/"{gene_name}.mafft.evalmsa.fasta"
    output:
        trimmed = PATHCURATED/"{gene_name}.mafft.evalmsa.trimmal.fasta"
    log:
        LOGDIR/"trimal"/"{gene_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        trimal -in {input.aligned} -out {wildcards.gene_name} -gt {params.threshold} 
        sed 's/ .*//g' {wildcards.gene_name} > {output.trimmed}
        """
    

