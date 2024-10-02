rule clean_headers:
    input:
        fasta = TARGET_DIR/("{gene_name}%s" % config["EXTENSION"])
    output:
        temp(PATHALN/"{gene_name}.clean.fasta")
    shell:
        """
        sed -e 's/_[0-9]\{{5\}}\($\| .*\)//' \
            -e 's/_R_//g' \
            -e 's/;.*//g' {input.fasta} > {output}
        """

rule align:
    threads: 8
    conda: "../envs/mafft.yaml"
    shadow: "minimal"
    input:
        fasta = PATHALN/"{gene_name}.clean.fasta"
    output:
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    shell:
        """
        mafft --thread {threads} \
            --adjustdirection {input.fasta} > {output.aligned}
        """

rule generate_consensus:
    conda:
        "../envs/biopython.yaml"
    input:
        aligned = PATHALN/"{gene_name}.mafft.fasta"
    output:
        consensus = PATHCONS/"{gene_name}.consensus.fasta"
    script:
        "../scripts/consensus.py"

rule add_consensus:
    input:
        consensus = PATHCONS/"{gene_name}.consensus.fasta",
        aln = PATHALN/"{gene_name}.mafft.fasta"
    output:
        with_consensus = PATHCONS/"{gene_name}.mafft.concat.fasta"
    shell:
        """
        sed -i 's/-/n/g' {input.consensus}
        cat {input.consensus} {input.aln} > {output.with_consensus}
        """

rule translate:
    threads: 8
    shadow: "shallow"
    conda:
        "../envs/fasta.yaml"
    input:
        PATHCONS/"{gene_name}.mafft.concat.fasta"
    output:
        translated = PATHCONS/"{gene_name}.mafft.concat.aa.fasta"
    shell:
        """
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
    input:
        translated = PATHCONS/"{gene_name}.mafft.concat.aa.fasta"
    output:
        vcf = PATHCONS/"{gene_name}.mafft.concat.aa.vcf"
    shell:
        """
        snp-sites -v {input.translated} > {output.vcf}
        """
