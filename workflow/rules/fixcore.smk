rule prune_trees:
    conda:
        "../envs/renv.yaml"
    params:
        max_missing = config["PRUNE_TREES"]["MAX_MISSING"],
        max_consecutve_variants = config["PRUNE_TREES"]["MAX_CONSECUTIVE_VARIANTS"]
    input:
        vcf = PATHCONS/"{gene_name}.mafft.concat.aa.vcf",
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    output:
        filtered = PATHCURATED/"{gene_name}.mafft.evalmsa.fasta",
        removed = PATHCONS/"{gene_name}.removed.names"
    log:
        LOGDIR/"prune_trees"/"{gene_name}.log"

    script:
        "../scripts/filter_alignments.R"
    
