function check_paths(
        data_dir::String;
        remake::Bool,
        n_job::Int=1,
    )

    println("Checking paths ...")

    dna_fa_gz = "$data_dir/grch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

    chromosome_bed_gz = "$data_dir/grch/chromosome.bed.gz"

    chrn_n_tsv = "$data_dir/grch/chrn_n.tsv"

    cdna_fa_gz = "$data_dir/grch/Homo_sapiens.GRCh38.cdna.all.fa.gz"

    virus_cdna_fa_gz = "$data_dir/virus/sequences.fasta"

    for path in (
            dna_fa_gz,
            chromosome_bed_gz,
            chrn_n_tsv,
            cdna_fa_gz,
            "$data_dir/grch/enst_gene_name.tsv",
            virus_cdna_fa_gz,
            "$data_dir/virus/sequences.csv",
        )

        if ispath(path)

            println(path)

        else

            error("$path doesn't exist.")

        end

    end

    dna_fa_bgz = replace(
        dna_fa_gz,
        ".gz"=>".bgz",
    )

    if remake || !ispath(dna_fa_bgz)

        Kraft.print_and_run_cmd(pipeline(
                `gzip --decompress $dna_fa_gz --stdout`,
                `bgzip --threads $n_job --stdout`,
                dna_fa_bgz,
                ))

    end

    dna_fa_gz_mmi = dna_fa_gz * ".mmi"

    if remake || !ispath(dna_fa_gz_mmi)

        Kraft.print_and_run_cmd(`minimap2 -t $n_job $dna_fa_gz -d $dna_fa_gz_mmi`)

    end

    if remake || !(ispath(dna_fa_bgz * ".fai") && ispath(dna_fa_bgz * ".gzi"))

        Kraft.print_and_run_cmd(`samtools faidx $dna_fa_bgz`)

    end

    if remake || !ispath(chromosome_bed_gz * ".tbi")

        Kraft.print_and_run_cmd(`tabix --force $chromosome_bed_gz`)

    end

    cdna_fa_gz_kallisto_index = cdna_fa_gz * ".kallisto_index"

    if remake || !ispath(cdna_fa_gz_kallisto_index)

        Kraft.print_and_run_cmd(`kallisto index $cdna_fa_gz --index $cdna_fa_gz_kallisto_index`)

    end

    virus_cdna_fa_gz_kallisto_index = virus_cdna_fa_gz * ".kallisto_index"

    if remake || !ispath(virus_cdna_fa_gz_kallisto_index)

        Kraft.print_and_run_cmd(`kallisto index $virus_cdna_fa_gz --index $virus_cdna_fa_gz_kallisto_index`)

    end

end
