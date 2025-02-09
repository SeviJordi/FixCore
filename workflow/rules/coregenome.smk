
rule panaroo:
    threads: config["CORE"]["N_CORES"]
    conda: "../envs/coregenome.yaml"
    params:
        threshold = config["CORE"]["THRESHOLD"],
        family_threshold = config["CORE"]["PANAROO"]["FAMILY_THRESHOLD"],
        clean_mode = config["CORE"]["PANAROO"]["CLEAN_MODE"],
        extra_params = config["CORE"]["PANAROO"]["EXTRA_PARAMS"]
    input:
        gff = expand(PROKKA_DIR/"{genome_name}/{genome_name}.gff", genome_name=iter_genome_names())
    output:
        out = directory(OUTDIR/config["CORE"]["TOOL"])
    shell:
        """
        panaroo -i {input.gff} \
            -o {output.out} \
            --clean-mode {params.clean_mode} \
            -t {threads} \
            --alignment core \
            --aligner mafft \
            --core_threshold {params.threshold} \
            -f {params.family_threshold} \
            {params.extra_params}
        """


rule select_core_genes:
    params:
        pantool = config["CORE"]["TOOL"],
        tar_dir = TARGET_DIR,
    input:
        indir = directory(OUTDIR/config["CORE"]["TOOL"])
    output:
        files = dynamic(TARGET_DIR/"{gene_name}.fasta")
    shell:
        """
        if [[{params.pantool} == "panacota"]]; then
            for file in {input.indir}/**/.gen; do
                name=$(basename $file .gen | sed 's/current.//')
                cp $file {params.tar_dir}/$name.fasta
        
        elif [[{params.pantool} == "panaroo"]]; then
            for file in {input.indir}/aligned_gene_sequences/*.aln.fas; do
                name=$(basename $file .aln.fas)
                cp $file {params.tar_dir}/$name.fasta

        elif [[{params.pantool} == "roary"]]; then
            gene_names=$(grep label {input.indir}/core_alignment_header.embl | cut -d"=" -f2)
            files=($(parallel echo {input.indir}/pan_genome_sequences/{{}}.fa.aln ::: $gene_names))
            for file in ${files[@]}; do 
                name=$(basename $file .fa.aln)
                cp $file {params.tar_dir}/$name.fasta
        """
