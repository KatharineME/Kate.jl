include("print_and_run_cmd.jl")


function make_path_dict(data_dir::String)

    println("Making path dict ...")

    dna_fa_gz::String = joinpath(
        data_dir,
        "grch",
        "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    )

    chromosome_bed_gz::String = joinpath(
        data_dir,
        "grch",
        "chromosome.bed.gz",
    )

    chrn_n_tsv::String = joinpath(
        data_dir,
        "grch",
        "chrn_n.tsv",
    )

    cdna_fa_gz::String = joinpath(
        data_dir,
        "grch",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
    )

    grch_enst_gene_name_tsv::String = joinpath(
        data_dir,
        "grch",
        "enst_gene_name.tsv",
    )

    virus_cdna_fa_gz::String = joinpath(
        data_dir,
        "virus",
        "sequences.fasta",
    )

    virus_sequences_csv::String = replace(
        virus_cdna_fa_gz,
        ".tsv"=>".csv",
    )

    for path::String in (
        dna_fa_gz,
        chromosome_bed_gz,
        chrn_n_tsv,
        cdna_fa_gz,
        grch_enst_gene_name_tsv,
        virus_cdna_fa_gz,
        virus_sequences_csv,
    )

        if ispath(path)

            println(path)

        else

            error("$path doesn't exist.")

        end

    end

    dna_fa_bgz::String = replace(
        dna_fa_gz,
        ".gz"=>".bgz",
    )

    if !ispath(dna_fa_bgz)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fa_gz --stdout`,
            `bgzip --threads 2 --stdout`,
            dna_fa_bgz,
        ))

    end

    Dict(
        "dna.fa.gz"=>dna_fa_gz,
        "dna.fa.bgz"=>dna_fa_bgz,
        "chromosome.bed.gz"=>chromosome_bed_gz,
        "chrn_n.tsv"=>chrn_n_tsv,
        "cdna.fa.gz"=>cdna_fa_gz,
        "virus_cdna.fa.gz"=>virus_cdna_fa_gz,
    )

end
