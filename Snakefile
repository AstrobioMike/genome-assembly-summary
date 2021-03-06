configfile: "config.yaml"

genome_IDs, = glob_wildcards(config["genomes_dir"] + "/{id}" + config["assembly_extension"])

# checking any were found
if len(genome_IDs) == 0:
    print("\n    No assembly fasta files were found, is the config.yaml file set appropriately?\n")
    exit(1)

# making directory for log files
try:
    os.mkdir(config["logs_dir"])
except:
    pass

rule all:
    input:
        str(config["output_prefix"]) + "-genome-summaries.tsv"
    shell:
        """
        rm -rf checkm-output/ checkm-summary.tsv gtdb-taxonomy.tsv gtdb-tk-out/ *-prots.faa
        """


rule gen_summary_stats:
    conda:
        "envs/bit.yaml"
    input:
        expand(config["genomes_dir"] + "/{genome_ID}"  + config["assembly_extension"], genome_ID = genome_IDs)
    output:
        "summary-stats.tsv"
    shell:
        """
        bit-summarize-assembly {input} -o {output}
        """


## processing if eukarya
if config["is_euk"]:

    rule run_eukcc:
        conda:
            "envs/eukcc.yaml"
        input:
            genome = config["genomes_dir"] + "/{genome_ID}" + config["assembly_extension"],
            eukcc_db_trigger = config["DIR_HOLDING_eukcc_DIR"] + "/" + config["eukcc_db_dir"] + "/" + config["eukcc_TRIGGER_FILE"]
        params:
            output_dir = "{genome_ID}-eukcc-out",
            eukcc_db_dir = config["DIR_HOLDING_eukcc_DIR"] + "/" + config["eukcc_db_dir"]
        resources:
            cpus = config["threads"]
        log:
            config["logs_dir"] + "{genome_ID}-eukcc.log"
        output:
            est_tab = "{genome_ID}-eukcc-estimates.tsv",
            AA_seqs = "{genome_ID}-prots.faa"
        shell:
            """
            eukcc single --db {params.eukcc_db_dir} --threads {resources.cpus} --out {params.output_dir} --keep {input.genome} > {log} 2>&1

            # getting comp./redund. estimates
            paste <( echo "{wildcards.genome_ID}" ) <( tail -n +2 {params.output_dir}/eukcc.csv | head -n 1 | cut -f 2,3 ) > {output.est_tab}

            # getting protein seqs to pass to CAT formatted in a way it likes so they can be tracked back to the contigs
            sed 's/metaeuk_//' {params.output_dir}/workdir/metaeuk/*_metaeuk_cleaned.faa > {output.AA_seqs}

            rm -rf {params.output_dir}
            """


    rule combine_eukcc_estimates:
        input:
            expand("{genome_ID}-eukcc-estimates.tsv", genome_ID = genome_IDs)
        output:
            "combined-eukcc-estimates.tsv"
        shell:
            """
            printf "Assembly\\tEst. Comp.\\tEst. Redund.\\n" > {output}
            cat {input} >> {output}
            rm {input}
            """


    rule run_CAT:
        conda:
            "envs/cat.yaml"
        input:
            genome = config["genomes_dir"] + "/{genome_ID}" + config["assembly_extension"],
            AA_seqs = "{genome_ID}-prots.faa",
            cat_db_trigger = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DIR"] + "/" + config["CAT_TRIGGER_FILE"]
        params:
            tmp_out_prefix = "{genome_ID}-tax-dir.tmp",
            tmp_tax = "{genome_ID}-tax.tmp",
            cat_db = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DIR"] + "/" + config["CAT_DB"],
            cat_tax = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DIR"] + "/" + config["CAT_TAX"],
            num_threads = config["threads"],
            assembly_extension = config["assembly_extension"]
        log:
            config["logs_dir"] + "{genome_ID}-CAT.log"
        output:
            "{genome_ID}-tax.tsv"
        shell:
            """
            CAT bin -b {input.genome} -p {input.AA_seqs} -d {params.cat_db} -t {params.cat_tax} -n {params.num_threads} -o {params.tmp_out_prefix} -f 0.5 -r 3 --top 4 --I_know_what_Im_doing > {log} 2>&1

            # adding names
            CAT add_names -i {params.tmp_out_prefix}.bin2classification.txt -o {params.tmp_tax} -t {params.cat_tax} --only_official > {log} 2>&1

            # formatting classification
            grep -v "^#" {params.tmp_tax} | awk -F $'\\t' ' BEGIN {{ OFS=FS }} {{ if ( $2 == "taxid assigned" ) {{ print $1,$6,$7,$8,$9,$10,$11,$12 }} \
            else {{ print $1,"NA","NA","NA","NA","NA","NA","NA" }} }} ' | head -n 1 | \
            sed 's/: [0-9\.]*//g' | sed 's/not classified/NA/g' | sed 's/no support/NA/g' | sed 's/{params.assembly_extension}//' > {output}

            rm -rf {wildcards.genome_ID}*tmp* {input.AA_seqs}
            """


    rule combine_euk_tax_outputs:
        input:
            expand("{genome_ID}-tax.tsv", genome_ID = genome_IDs)
        output:
            "CAT-taxonomies.tsv"
        shell:
            """
            printf "Assembly\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n" > {output}
            cat {input} >> {output}
            rm {input}
            """


    rule combine_euk_outputs:
        input:
            assembly_stats_tab = "summary-stats.tsv",
            eukcc_tab = "combined-eukcc-estimates.tsv",
            tax_tab = "CAT-taxonomies.tsv"
        output:
            str(config["output_prefix"]) + "-genome-summaries.tsv"
        shell:
            """
            python scripts/combine-euk-outputs.py -s {input.assembly_stats_tab} -c {input.eukcc_tab} -t {input.tax_tab} -o {output}
            rm {input}
            """

# processing if not eukarya
if not config["is_euk"]:
    rule run_checkm:
        conda:
            "envs/checkm.yaml"
        params:
            pplacer_cpus = config["gtdb_tk_checkm_pplacer_cpus"],
            genomes_dir = config["genomes_dir"],
            extension = config["assembly_extension"]
        resources:
            cpus = config["threads"]
        output:
            "checkm-summary.tsv"
        log:
            config["logs_dir"] + "checkm.log"
        shell:
            """
            checkm lineage_wf -x {params.extension} -t {resources.cpus} --pplacer_threads {params.pplacer_cpus} --tab_table -f checkm-summary.tsv {params.genomes_dir} checkm-output > {log} 2>&1
            """


    rule gtdb_tk_classify:
        conda:
            "envs/gtdb-tk.yaml"
        input:
            gtdbtk_db_trigger = config["GTDB_DATA_PATH"] + "/" + config["GTDB_TRIGGER_FILE"]
        params:
            pplacer_cpus = config["gtdb_tk_checkm_pplacer_cpus"],
            genomes_dir = config["genomes_dir"],
            extension = config["assembly_extension"]
        resources:
            cpus = config["threads"]
        output:
            directory("gtdb-tk-out")
        log:
            config["logs_dir"] + "gtdb.log"
        shell:
            """
            gtdbtk classify_wf -x {params.extension} --genome_dir {params.genomes_dir} --out_dir {output} --cpus {resources.cpus} --pplacer_cpus {params.pplacer_cpus} > {log} 2>&1
            """


    rule combine_gtdb_classifications:
        input:
            "gtdb-tk-out"
        output:
            "gtdb-taxonomy.tsv"
        shell:
            "cut -f 1,2 gtdb-tk-out/classify/*summary.tsv > {output}"


    rule combine_outputs:
        input:
            assembly_stats_tab = "summary-stats.tsv",
            checkm_tab = "checkm-summary.tsv",
            gtdb_tax_tab = "gtdb-taxonomy.tsv"
        output:
            str(config["output_prefix"]) + "-genome-summaries.tsv"
        shell:
            "python scripts/combine-outputs.py -s {input.assembly_stats_tab} -c {input.checkm_tab} -t {input.gtdb_tax_tab} -o {output}"

rule clean:
    shell:
        "rm -rf checkm-output checkm-snake.log checkm-summary.tsv gtdb-snake.log gtdb-taxonomy.tsv gtdb-tk-out summary-stats.tsv *-genome-summaries.tsv combined-eukcc-estimates.tsv CAT-taxonomies.tsv *-renamed-prots.faa *-CAT.log *-eukcc.log"


### database checking and setup rules ###

rule setup_gtdbtk_db:
    """
    This rule checks for the gtdb-tk db (minimally currently) and downloads if needed.
    """

    conda:
        "envs/gtdb-tk.yaml"
    output:
        gtdbtk_db_trigger = config["GTDB_DATA_PATH"] + "/" + config["GTDB_TRIGGER_FILE"]
    params:
        gtdbtk_db_dir = config["GTDB_DATA_PATH"]
    log:
        config["logs_dir"] + "setup-gtdbtk-db.log"
    shell:
        """
        mkdir -p {params.gtdbtk_db_dir}
        # storing current working directory to be able to send the log file here
        working_dir=$(pwd)
        cd {params.gtdbtk_db_dir}
        # adding wanted location to this conda env PATH (gtdb-tk looks in the GTDBTK_DATA_PATH variable),
            # so will be set when the conda environment is started from now on
        mkdir -p ${{CONDA_PREFIX}}/etc/conda/activate.d/
        echo 'export GTDBTK_DATA_PATH={params.gtdbtk_db_dir}' >> ${{CONDA_PREFIX}}/etc/conda/activate.d/set_env_vars.sh
        # but still needs to be set for this particular session that is downloading and setting up the db
        GTDBTK_DATA_PATH={params.gtdbtk_db_dir}
        # now downloading
        download-db.sh > ${{working_dir}}/{log} 2>&1
        cd - > /dev/null
        touch {output.gtdbtk_db_trigger}
        """


rule setup_CAT_db:
    """
    This rule checks for the CAT reference database, and downloads if needed.
    """

    conda:
        "envs/cat.yaml"
    output:
        cat_db_trigger = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DIR"] + "/" + config["CAT_TRIGGER_FILE"]
    params:
        cat_db_dir = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DIR"],
        compressed_cat = config["DIR_HOLDING_CAT_DIR"] + "/" + config["CAT_DL_FILE"],
        compressed_nr_faa = config["CAT_DIR"] + "/" + config["CAT_DB"] + "/" + config["CAT_COMPRESSED_NR_FAA"],
        cat_dl_link = config["CAT_DL_LINK"],
        DIR_HOLDING_CAT_DIR = config["DIR_HOLDING_CAT_DIR"]
    log:
        config["logs_dir"] + "setup-CAT-db.log"
    shell:
        """
        mkdir -p {params.DIR_HOLDING_CAT_DIR}
        printf "### Setting up CAT reference database ###\n\n" > {log} 2>&1
        printf "  Downloading reference db:\n\n" >> {log} 2>&1
        curl -L -C - -o {params.compressed_cat} {params.cat_dl_link} >> {log} 2>&1
        printf "\n\n  Extracting reference db:\n\n" >> {log} 2>&1
        tar -xvzf {params.compressed_cat} -C {params.DIR_HOLDING_CAT_DIR} >> {log} 2>&1
        rm {params.compressed_cat} {params.compressed_nr_faa}
        touch {output.cat_db_trigger}
        """


rule setup_eukcc_db:
    """
    This rule checks for the eukcc database, and downloads if needed.
    """

    conda:
        "envs/eukcc.yaml"
    output:
        eukcc_db_trigger = config["DIR_HOLDING_eukcc_DIR"] + "/" + config["eukcc_db_dir"] + "/" + config["eukcc_TRIGGER_FILE"]
    params:
        compressed_eukcc = config["DIR_HOLDING_eukcc_DIR"] + "/" + config["eukcc_DL_FILE"],
        eukcc_dl_link = config["eukcc_DL_LINK"],
        DIR_HOLDING_eukcc_DIR = config["DIR_HOLDING_eukcc_DIR"]
    log:
        config["logs_dir"] + "setup-eukcc-db.log"
    shell:
        """
        mkdir -p {params.DIR_HOLDING_eukcc_DIR}
        printf "### Setting up eukcc reference database ###\n\n" > {log} 2>&1
        printf "  Downloading reference db:\n\n" >> {log} 2>&1
        curl -L -C - -o {params.compressed_eukcc} {params.eukcc_dl_link} >> {log} 2>&1
        printf "\n\n  Extracting reference db:\n\n" >> {log} 2>&1
        tar -xvzf {params.compressed_eukcc} -C {params.DIR_HOLDING_eukcc_DIR} >> {log} 2>&1
        rm {params.compressed_eukcc}
        touch {output.eukcc_db_trigger}
        """
