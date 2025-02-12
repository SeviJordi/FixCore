
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
        flag = OUTDIR/".panaroo_done"
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
        family_threshold = config["CORE"]["ROARY"]["FAMILY_THRESHOLD"]*100 if config["CORE"]["ROARY"]["FAMILY_THRESHOLD"] < 1 else config["CORE"]["ROARY"]["FAMILY_THRESHOLD"],
        extra_params = config["CORE"]["ROARY"]["EXTRA_PARAMS"],
        out = PATHPAN
    input:
        gff = expand(PROKKA_DIR/"{genome_name}/{genome_name}.gff", genome_name=iter_genome_names())
    output:  
        flag = OUTDIR/".roary_done"
    log:
        LOGDIR/"roary"/"roary.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        touch {output.flag}
        
        [[ -d {params.out} ]] && rm -r {params.out}
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
        outdir = PATHPAN,
        name = config["CORE"]["PANACOTA"]["NAME"]
    input:
        fastas = expand(GENOMES_DIR/"{genome_name}.fasta", genome_name=iter_genome_names())
    output:
        PATHPAN/"sample_sheet.lst"
    log:
        LOGDIR/"panacota"/"sample_sheet.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        mkdir -p {params.outdir}

        for file in {input.fastas}; do
            echo "$( basename $file) :: {params.name}"
        done > {output}
        """

rule panacota_annotate:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        out = PATHPAN,
        genomes = GENOMES_DIR,
        name = config["CORE"]["PANACOTA"]["NAME"],
        nbcont = config["CORE"]["PANACOTA"]["ANNOTATE"]["NBCONT"],
        l90 = config["CORE"]["PANACOTA"]["ANNOTATE"]["L90"],
        cutn = config["CORE"]["PANACOTA"]["ANNOTATE"]["CUTN"],
        extra_params = config["CORE"]["PANACOTA"]["ANNOTATE"]["EXTRA_PARAMS"]
    input:
        sample_sheet = PATHPAN/"sample_sheet.lst"
    output:
        listfile = PATHPAN/"annotation/LSTINFO-sample_sheet.lst",
    log:
        LOGDIR/"panacota"/"annotate.log"
    shell:
        """
        exec >{log}
        exec 2>&1
        pwd
        run_panacota.py annotate --nbcont {params.nbcont} \
            --l90 {params.l90} \
            --cutn {params.cutn} \
            -n {params.name} \
            -r {params.out}/annotation \
            --threads {threads} \
            -l {input.sample_sheet} \
            -q \
            -d {params.genomes} \
            {params.extra_params}
        """  

rule panacota_pangenome:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        out = PATHPAN,
        prefix = PREFIX,
        extra_params = config["CORE"]["PANACOTA"]["PANGENOME"]["EXTRA_PARAMS"],
        family_threshold = config["CORE"]["PANACOTA"]["FAMILY_THRESHOLD"]
    input:
        listfile = PATHPAN/"annotation/LSTINFO-sample_sheet.lst"
    output:
        pangenome = PATHPAN/f"pangenome/{PREFIX}_pangenome.lst",
    log:
        LOGDIR/"panacota"/"pangenome.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        run_panacota.py pangenome -l {input.listfile} \
            -n {params.prefix} \
            -d {params.out}/annotation/Proteins \
            -o {params.out}/pangenome \
            --threads {threads} \
            -i {params.family_threshold} \
            -q \
            -f $( basename {output.pangenome}) \
            {params.extra_params}
        """

rule panacota_corepers:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        threshold = config["CORE"]["THRESHOLD"],
        extra_params = config["CORE"]["PANACOTA"]["COREPERS"]["EXTRA_PARAMS"],
        mode = config["CORE"]["PANACOTA"]["COREPERS"]["MODE"]
    input:
        pangenome = PATHPAN/f"pangenome/{PREFIX}_pangenome.lst"
    output:
        dir_corepers = directory(PATHPAN/f"corepers/")
    log:
        LOGDIR/"panacota"/"corepers.log"
    shell:
        """
        exec >{log}
        exec 2>&1
        
        if [[ "{params.mode}" == "mixed" ]]; then
            flag="-X"
        elif [[ "{params.mode}" == "multi" ]]; then
            flag="-M"
        elif [[ "{params.mode}" == "" ]]; then
            flag=""
        else
            echo "ERROR: Select a valid mode for core genome constriction: mixed or multi"
            exit 1
        fi

        run_panacota.py corepers -p {input.pangenome} \
            -t {params.threshold} \
            -o {output.dir_corepers} \
            $flag {params.extra_params}

        """

rule panacota_align:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/panacota.yaml"
    params:
        out = PATHPAN,
        prefix = PREFIX,
        threshold = config["CORE"]["THRESHOLD"],
        extra_params = config["CORE"]["PANACOTA"]["ALIGN"]["EXTRA_PARAMS"]
    input:
        dir_corepers = PATHPAN/f"corepers/",
        listfile = PATHPAN/"annotation/LSTINFO-sample_sheet.lst"
    output:
        directory(PATHPAN/f"alignment/Align-{PREFIX}")
    log:
        LOGDIR/"panacota"/"align.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        run_panacota.py align --threads {threads} \
            -F \
            -n {params.prefix} \
            -c {input.dir_corepers}/PersGenome*{params.threshold}*.lst \
            -l {input.listfile} \
            -d {params.out}/annotation \
            -o {params.out}/alignment \
            {params.extra_params}
        """

rule rename_panacota:
    conda: "../envs/biopython.yaml"
    input:
        dir_al = PATHPAN/f"alignment/Align-{PREFIX}",
        listfile = PATHPAN/"annotation/LSTINFO-sample_sheet.lst"
    output:
        done = OUTDIR/".panacota_done"
    script:
        "../scripts/rename_panacota.py"
        
checkpoint select_core_genes:
    conda: "../envs/common.yaml"
    params:
        pantool = config["CORE"]["TOOL"],
        tar_dir = TARGET_DIR,
        indir = OUTDIR/config["CORE"]["TOOL"],
        fix_genes = PATHCURATED,
        concatenated = OUTDIR/f"{PREFIX}.concatenated.fixcore.fasta"
    input:
        done = OUTDIR/f".{config['CORE']['TOOL']}_done"
    output:
        directory(TARGET_DIR)
    log:
        LOGDIR/"select_core_genes"/"select_core_genes.log"
    shell:
        """
        exec >{log}
        exec 2>&1

        mkdir -p {params.tar_dir}
        rm {input.done}

        if ls {params.tar_dir}*.fasta &> /dev/null; then
            rm {params.tar_dir}/*
        fi

        if ls {params.fix_genes}*.fasta &> /dev/null; then
            rm {params.fix_genes}/*
        fi

        [[ -f {params.concatenated} ]] && rm {params.concatenated}

        if [[ "{params.pantool}" == "panacota" ]]; then
            for file in {params.indir}/alignment/**/*.gen; do
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
