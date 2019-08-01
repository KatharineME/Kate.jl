using Dates

include("print_and_run_cmd.jl")


function trim_sequence(
    _1_fastq_gz_file_path::String,
    _2_fastq_gz_file_path::String,
    output_file_path_prefix::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Trimming sequence ...")

    if isfile("$output_file_path_prefix-trimmed-pair1.fastq.gz") || isfile("$output_file_path_prefix-trimmed-pair2.fastq.gz")

        error("$output_file_path_prefix-trimmed-pair(1|2).fastq.gz exists.")

    end
    
    output_directory_path = splitdir(output_file_path_prefix)[1]

    mkpath(output_directory_path)

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_file_path_prefix --quiet $_1_fastq_gz_file_path $_2_fastq_gz_file_path`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
