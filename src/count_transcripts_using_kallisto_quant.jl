include("print_command_and_run.jl")

function count_transcripts_using_kallisto_quant(
    fastq_gz_file_paths::Array,
    fasta_gz_file_path::String,
    output_directory_path::String;
    fragment_length::Int=180,
    fragment_length_standard_deviation::Int=20,
    n_job::Int=1,
    overwrite::Bool=false,
)::String

    fasta_gz_kallisto_index_file_path = "$fasta_gz_file_path.kallisto.index"

    if !isfile(fasta_gz_kallisto_index_file_path)

        print_command_and_run(`kallisto index --index $fasta_gz_kallisto_index_file_path $fasta_gz_file_path`)

    end

    abundance_file_path = "$output_directory_path/abundance.tsv"

    if !overwrite && isfile(abundance_file_path)

        error(abundance_file_path)

    end

    if length(fastq_gz_file_paths) == 1

        arguments = [
            "--single",
            "--fragment-length",
            fragment_length,
            "--sd",
            fragment_length_standard_deviation,
        ]

    elseif length(fastq_gz_file_paths) == 2

        arguments = []

    end

    print_command_and_run(`kallisto quant --index $fasta_gz_kallisto_index_file_path --output-dir $output_directory_path --threads $n_job $arguments $fastq_gz_file_paths`)

    return output_directory_path

end
