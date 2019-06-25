include("print_and_run_cmd.jl")


function check_data_for_processing_sequence(
    data_dir::String,
    n_job::Int,
)

    println("Check data for processing sequence ...")

    dna_fasta_gz::String = joinpath(
        data_dir,
        "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    )

    chromosome_bed_gz::String = joinpath(
        data_dir,
        "chromosome.bed.gz",
    )

    chrn_n_tsv::String = joinpath(
        data_dir,
        "chrn_n.tsv",
    )

    cdna_fasta_gz::String = joinpath(
        data_dir,
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
    )

    enst_gene_name_tsv::String = joinpath(
        data_dir,
        "enst_gene_name.tsv",
    )

    for path::String in (
        dna_fasta_gz,
        chromosome_bed_gz,
        chrn_n_tsv,
        cdna_fasta_gz,
        enst_gene_name_tsv,
    )

        if isfile(path)

            println(path)

        else

            error("$path doesn't exist.")

        end

    end

    dna_fasta_bgz::String = "$(splitext(dna_fasta_gz)[1]).bgz"

    if !isfile(dna_fasta_bgz)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fasta_gz --stdout`,
            `bgzip --threads $n_job --stdout`,
            dna_fasta_bgz,
        ))

    end

end
