using Dates

include("print_and_run_cmd.jl")


function trim_sequence(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    output_prefix::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Trimming sequence ...")

    output_dir = splitdir(output_prefix)[1]

    if isdir(output_dir)

        error("$output_dir exists.")
    
    else

        mkpath(output_dir)
    
    end

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end
