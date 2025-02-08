
rule align:
    threads: config["MAFFT"]["N_CORES"]
    conda: "../envs/mafft.yaml"
    shadow: "minimal"
    input:
        fasta = TARGET_DIR/"{gene_name}.fasta"
    output:
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    log:
        LOGDIR/"align"/"{gene_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        mafft --thread {threads} \
            --adjustdirection {input.fasta} \
         sed -e 's/_[0-9]\{{5\}}\($\| .*\)//' \
            -e 's/_R_//g' \
            -e 's/;.*//g > {output.aligned}
        """

rule generate_consensus:
    conda:
        "../envs/biopython.yaml"
    threads: 1
    input:
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    output:
        consensus = PATHCONS/"{gene_name}.consensus.fasta"
    log:
        LOGDIR/"consensus"/"{gene_name}.log"
    script:
        "../scripts/consensus.py"

rule add_consensus:
    threads: 1
    input:
        consensus = PATHCONS/"{gene_name}.consensus.fasta",
        aln = PATHALN/"{gene_name}.mafft.fasta"
    output:
        with_consensus = PATHCONS/"{gene_name}.mafft.concat.fasta",
    log:
        LOGDIR/"add_consensus"/"{gene_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        sed -i 's/-/n/g' {input.consensus}
        cat {input.consensus} {input.aln} > {output.with_consensus}
        """

rule translate:
    threads: config["SEQKIT"]["N_CORES"]
    shadow: "shallow"
    conda:
        "../envs/fasta.yaml"
    input:
        PATHCONS/"{gene_name}.mafft.concat.fasta"
    output:
        translated = PATHCONS/"{gene_name}.mafft.concat.aa.fasta"
    log:
        LOGDIR/"translate"/"{gene_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        seqkit translate --threads {threads} {input} > tmp1
        seqkit replace --threads {threads} \
            -p "(-)" -r 'X' \
            -s tmp1 > tmp2
        seqkit replace --threads {threads} \
            -p "(N)" -r '*' \
            -s tmp2 > {output.translated}
        """

rule generate_vcf:
    conda:
        "../envs/fasta.yaml"
    threads: 1
    input:
        translated = PATHCONS/"{gene_name}.mafft.concat.aa.fasta"
    output:
        vcf = PATHCONS/"{gene_name}.mafft.concat.aa.vcf"
    log:
        LOGDIR/"generate_vcf"/"{gene_name}.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        snp-sites -v {input.translated} | cat > {output.vcf}
        """
