rule prune_trees:
    conda:
        "../envs/renv.yaml"
    input:
        vcf = PATHCONS/"{gene_name}.mafft.concat.aa.vcf",
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    output:
        filtered = PATHCURATED/"{gene_name}.mafft.evalmsa.fasta",
        removed = PATHCONS/"{gene_name}.removed.names"

    script:
        "../scripts/filter_alignments.R"
    
