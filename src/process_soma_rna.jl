include("trim_sequence.jl")

include("check_sequence.jl")

include("count_transcript.jl")


function process_soma_rna(
    soma_rna_1_fastq_gz::String,
    soma_rna_2_fastq_gz::String,
    output_dir::String,
    cdna_fasta_gz::String,
    n_job::Int,
)

    for file_path in(
        soma_rna_1_fastq_gz,
        soma_rna_2_fastq_gz,
        cdna_fasta_gz,
    )

        if !isfile(file_path)

            error("$file_path doesn't exist.")

        end

    end

    soma_trim_sequence_prefix = joinpath(
        output_dir,
        "trim_sequence",
        "soma",
    )

    trim_sequence(
        soma_rna_1_fastq_gz,
        soma_rna_2_fastq_gz,
        soma_trim_sequence_prefix,
        n_job,
    )

    soma_trim_1_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    soma_trim_2_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    check_sequence(
        (
            soma_trim_1_fastq_gz,
            soma_trim_2_fastq_gz,
        ),
        joinpath(
            output_dir,
            "check_sequence",
        ),
        n_job,
    )

    count_transcript_dir = joinpath(
        output_dir,
        "count_transcript",
    )

    count_transcript(
        soma_trim_1_fastq_gz,
        soma_trim_2_fastq_gz,
        cdna_fasta_gz,
        count_transcript_dir,
        n_job,
    )

end
