include("align_sequence.jl")
include("check_sequence.jl")
include("find_variant.jl")
include("print_and_run_cmd.jl")
include("trim_sequence.jl")


function process_soma_dna(
    germ_dna_1_fastq_gz_file_path::String,
    germ_dna_2_fastq_gz_file_path::String,
    soma_dna_1_fastq_gz_file_path::String,
    soma_dna_2_fastq_gz_file_path::String,
    is_targeted::Bool,
    output_directory_path::String,
    dna_fasta_gz_file_path::String,
    chromosome_bed_gz_file_path::String,
    chrn_n_tsv_file_path::String,
    n_job::Int,
    n_gb_memory::Int,
    n_gb_memory_per_job::Int,
)

    for file_path in (
        germ_dna_1_fastq_gz_file_path,
        germ_dna_2_fastq_gz_file_path,
        soma_dna_1_fastq_gz_file_path,
        soma_dna_2_fastq_gz_file_path,
        dna_fasta_gz_file_path,
        chromosome_bed_gz_file_path,
        chrn_n_tsv_file_path,
    )

        if !isfile(file_path)

            error("$file_path does not exist.")

        end

    end

    if isdir(output_directory_path)

        error("$output_directory_path exists.")

    end

    mkpath(output_directory_path)

    germ_trim_sequence_file_path_prefix = joinpath(
        output_directory_path,
        "trim_sequence",
        "germ",
    )

    trim_sequence(
        germ_dna_1_fastq_gz_file_path,
        germ_dna_2_fastq_gz_file_path,
        germ_trim_sequence_file_path_prefix,
        n_job,
    )

    germ_trim_1_fastq_gz_file_path = "$germ_trim_sequence_file_path_prefix-trimmed-pair1.fastq.gz"

    germ_trim_2_fastq_gz_file_path = "$germ_trim_sequence_file_path_prefix-trimmed-pair2.fastq.gz"

    soma_trim_sequence_file_path_prefix = joinpath(
        output_directory_path,
        "trim_sequence",
        "soma",
    )

    trim_sequence(
        soma_dna_1_fastq_gz_file_path,
        soma_dna_2_fastq_gz_file_path,
        soma_trim_sequence_file_path_prefix,
        n_job,
    )

    soma_trim_1_fastq_gz_file_path = "$soma_trim_sequence_file_path_prefix-trimmed-pair1.fastq.gz"

    soma_trim_2_fastq_gz_file_path = "$soma_trim_sequence_file_path_prefix-trimmed-pair2.fastq.gz"

    check_sequence(
        (
         germ_trim_1_fastq_gz_file_path,
         germ_trim_2_fastq_gz_file_path,
         soma_trim_1_fastq_gz_file_path,
         soma_trim_2_fastq_gz_file_path,
        ),
        joinpath(output_directory_path, "check_sequence"),
        n_job,
    )

    germ_bam_file_path = joinpath(output_directory_path, "align_sequence", "germ.bam")

    align_sequence(
        germ_trim_1_fastq_gz_file_path,
        germ_trim_2_fastq_gz_file_path,
        "Germ",
        dna_fasta_gz_file_path,
        germ_bam_file_path,
        n_job,
        n_gb_memory_per_job,
    )

    soma_bam_file_path = joinpath(output_directory_path, "align_sequence", "soma.bam")

    align_sequence(
        soma_trim_1_fastq_gz_file_path,
        soma_trim_2_fastq_gz_file_path,
        "Soma",
        dna_fasta_gz_file_path,
        soma_bam_file_path,
        n_job,
        n_gb_memory_per_job,
    )

    dna_fasta_file_path_prefix = splitext(dna_fasta_gz_file_path)[1]

    dna_fasta_bgz_file_path = "$dna_fasta_file_path_prefix.bgz"

    if !isfile(dna_fasta_bgz_file_path)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fasta_gz_file_path --stdout`,
            `bgzip --threads $n_job --stdout`,
            dna_fasta_bgz_file_path,
        ))

    end

    find_variant_directory_path = joinpath(output_directory_path, "find_variant")

    find_variant(
        germ_bam_file_path,
        soma_bam_file_path,
        is_targeted,
        dna_fasta_bgz_file_path,
        chromosome_bed_gz_file_path,
        chrn_n_tsv_file_path,
        find_variant_directory_path,
        n_job,
        n_gb_memory,
    )

    return nothing

end
