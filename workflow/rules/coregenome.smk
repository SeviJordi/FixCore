
rule panaroo:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panaroo.yaml"
    params:
        threshold = config["CORE"]["THRESHOLD"],
        family_threshold = config["CORE"]["PANAROO"]["FAMILY_THRESHOLD"],
        clean_mode = config["CORE"]["PANAROO"]["CLEAN_MODE"],
        extra_params = config["CORE"]["PANAROO"]["EXTRA_PARAMS"],
        out = PATHPAN
    input:
        gff = expand(PROKKA_DIR/"{genome_name}/{genome_name}.gff", genome_name=iter_genome_names())
    output:  
        flag = temp(".panaroo_done")
    log:
        LOGDIR/"panaroo"/"panaroo.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        touch {output.flag}
        
        panaroo -i {input.gff} \
            -o {params.out} \
            --clean-mode {params.clean_mode} \
            -t {threads} \
            --alignment core \
            --aligner mafft \
            --core_threshold {params.threshold} \
            -f {params.family_threshold} \
            {params.extra_params}
        """

rule roary:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/roary.yaml"
    params:
        threshold = config["CORE"]["THRESHOLD"]*100 if config["CORE"]["THRESHOLD"] < 1 else config["CORE"]["THRESHOLD"],
        family_threshold = config["CORE"]["PANAROO"]["FAMILY_THRESHOLD"]*100 if config["CORE"]["PANAROO"]["FAMILY_THRESHOLD"] < 1 else config["CORE"]["PANAROO"]["FAMILY_THRESHOLD"],
        extra_params = config["CORE"]["ROARY"]["EXTRA_PARAMS"],
        out = PATHPAN
    input:
        gff = expand(PROKKA_DIR/"{genome_name}/{genome_name}.gff", genome_name=iter_genome_names())
    output:  
        flag = temp(".roary_done")
    log:
        LOGDIR/"roary"/"roary.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        touch {output.flag}
        
        roary -e \
            --mafft \
            -p {threads} \
            -i {params.family_threshold} \
            -f {params.out} \
            -v -z \
            -cd {params.threshold} \
            {params.extra_params} {input.gff}
        """


rule panacota_sample_sheet:
    params:
        outdir = PATHPAN
    input:
        fastas = expand(GENOMES_DIR/"{genome_name}.fasta", genome_name=iter_genome_names())
    output:
        PATHPAN/"sample_sheet.lst
    log:
        LOGDIR/"panacota"/"sample_sheet.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        mkdir -p {params.outdir}

        for file in {input.fastas}; do
            echo "$file :: $(basename $file .fasta)"
        done > {output}
        """

rule panacota_annotate:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        out = PATHPAN,
        genomes = GENOMES_DIR
    input:
        sample_sheet = PATHPAN/"sample_sheet.lst
    output:
        listfile = PATHPAN/"annotation/LSTINFO_sample_sheet.lst",
    log:
        LOGDIR/"panacota"/"annotate.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        run_panacota.py annotate --nbcont 2000 \
            --l90 500 \
            --cutn 0 \
            -r {params.out}/annotation \
            --threads {threads} \
            -l {input.sample_sheet} \
            -q \
            -d {params.genomes}
        """  

rule panacota_pangenome:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        out = PATHPAN,
        genomes = GENOMES_DIR
    input:
        listfile = PATHPAN/"annotation/LSTINFO_sample_sheet.lst"
    output:
    log:
        LOGDIR/"panacota"/"pangenome.log"
    shell:
        """
        exec >{log}
        exec 2>&1
        
        """

checkpoint select_core_genes:
    params:
        pantool = config["CORE"]["TOOL"],
        tar_dir = TARGET_DIR,
        indir = OUTDIR/config["CORE"]["TOOL"]
    input:
        done = f".{config['CORE']['TOOL']}_done"
    output:
        directory(TARGET_DIR)
    log:
        LOGDIR/"select_core_genes"/"select_core_genes.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        mkdir -p {params.tar_dir}

        if [[ "{params.pantool}" == "panacota" ]]; then
            for file in {params.indir}/**/.gen; do
                name=$(basename $file .gen | sed 's/current.//')
                cp $file {params.tar_dir}/$name.fasta
            done

        elif [[ "{params.pantool}"  == "panaroo" ]]; then
            for file in {params.indir}/aligned_gene_sequences/*.aln.fas; do
                name=$(basename $file .aln.fas)
                cp $file {params.tar_dir}/$name.fasta
            done

        elif [[ "{params.pantool}" == "roary" ]]; then
            gene_names=$(grep label {params.indir}/core_alignment_header.embl | cut -d"=" -f2)
            files=($(parallel echo {params.indir}/pan_genome_sequences/{{}}.fa.aln ::: $gene_names))
            for file in ${{files[@]}}; do 
                name=$(basename $file .fa.aln)
                cp $file {params.tar_dir}/$name.fasta
            done

        fi

        """
