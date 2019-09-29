include("check_sequence.jl")
include("count_transcript.jl")
include("trim_sequence.jl")


function process_soma_rna(
    soma_rna_1_fastq_gz_file_path::String,
    soma_rna_2_fastq_gz_file_path::String,
    output_directory_path::String,
    cdna_fasta_gz_file_path::String,
    n_job::Int,
)

    for file_path in (
        soma_rna_1_fastq_gz_file_path,
        soma_rna_2_fastq_gz_file_path,
        cdna_fasta_gz_file_path,
    )

        if !isfile(file_path)

            error("$file_path does not exist.")

        end

    end

    if isdir(output_directory_path)

        error("$output_directory_path exists.")

    end

    mkpath(output_directory_path)

    soma_trim_sequence_file_path_prefix = joinpath(
        output_directory_path,
        "trim_sequence",
        "soma",
    )

    trim_sequence(
        soma_rna_1_fastq_gz_file_path,
        soma_rna_2_fastq_gz_file_path,
        soma_trim_sequence_file_path_prefix,
        n_job,
    )

    soma_trim_1_fastq_gz_file_path = "$soma_trim_sequence_file_path_prefix-trimmed-pair1.fastq.gz"

    soma_trim_2_fastq_gz_file_path = "$soma_trim_sequence_file_path_prefix-trimmed-pair2.fastq.gz"

    check_sequence(
        (soma_trim_1_fastq_gz_file_path, soma_trim_2_fastq_gz_file_path),
        joinpath(output_directory_path, "check_sequence"),
        n_job,
    )

    count_transcript_directory_path = joinpath(output_directory_path, "count_transcript")

    count_transcript(
        soma_trim_1_fastq_gz_file_path,
        soma_trim_2_fastq_gz_file_path,
        cdna_fasta_gz_file_path,
        count_transcript_directory_path,
        n_job,
    )

    return nothing

end
