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

    if isfile("$output_prefix-trimmed-pair1.fastq.gz") || isfile("$output_prefix-trimmed-pair2.fastq.gz")

        error("$output_prefix-trimmed-pair(1|2).fastq.gz exists.")

    end
    
    output_dir = splitdir(output_prefix)[1]

    mkpath(output_dir)

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
