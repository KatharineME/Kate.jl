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

    trim_sequence_dir = splitdir(output_prefix)[1]

    germ_or_soma = splitdir(output_prefix)[2]

    trimmed_pair_1_path = joinpath(string(trim_sequence_dir), string(germ_or_soma, "-trimmed-pair1.fastq.gz"))

    trimmed_pair_2_path = joinpath(string(trim_sequence_dir), string(germ_or_soma, "-trimmed-pair2.fastq.gz"))

    if isfile(trimmed_pair_1_path) && isfile(trimmed_pair_2_path)
        
        println("Skipping trimming because trimmed files already exist:\n $trimmed_pair_1_path\n $trimmed_pair_2_path\n")

    else
        
        println("($start_time) Trimming sequence ...")
        
        output_dir::String = splitdir(output_prefix)[1]
        
        mkpath(output_dir)
        
        print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --mode pe -Q 2 -q 2 --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)
        
        # print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --end-quality 20 --mode pe --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)
        
        end_time = now()
        
        println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")        

    end

end
