using Dates

include("print_and_run_cmd.jl")


function count_transcript(
    _1_fastq_gz_file_path::String,
    _2_fastq_gz_file_path::String,
    cdna_fasta_gz_file_path::String,
    output_directory_path::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Counting transcript...")

    cdna_fasta_gz_kallisto_index_file_path = "$cdna_fasta_gz_file_path.kallisto_index"

    if !isfile(cdna_fasta_gz_kallisto_index_file_path)

        print_and_run_cmd(`kallisto index --index $cdna_fasta_gz_kallisto_index_file_path $cdna_fasta_gz_file_path`)

    end

    if isdir(output_directory_path)

        error("$output_directory_path exists.")

    end

    mkpath(output_directory_path)

    print_and_run_cmd(`kallisto quant --threads $n_job --index $cdna_fasta_gz_kallisto_index_file_path --output-dir $output_directory_path $_1_fastq_gz_file_path $_2_fastq_gz_file_path`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

    return nothing

end
