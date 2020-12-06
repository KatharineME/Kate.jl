using Dates

include("print_and_run_cmd.jl")


function trim_sequence(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    output_dir,
    output_prefix::String,
    n_job::Int,
)

    start_time = now()

    if isfile(joinpath(output_dir, "trim_sequence/germ-trimmed-pair1.fastq.gz")) && isfile(joinpath(output_dir, "trim_sequence/germ-trimmed-pair2.fastq.gz"))
        
        println("Trimmed files already exist.\n")

    else
        
        println("($start_time) Trimming sequence ...")
        
        output_dir::String = splitdir(output_prefix)[1]
        
        mkpath(output_dir)
        
        print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --end-quality 20 --mode pe --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)
        
        endtime = now()
        
        println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")        

    end

end
