using Dates

include("print_and_run_cmd.jl")


function check_sequence(fastq_gzs::Tuple{Vararg{String}}, output_dir::String, n_job::Int)

    start_time = now()

    if ispath(output_dir)

        println("Skipping check_sequence because $output_dir already exists.\n")

    else

        println("($start_time) Checking sequence ...")

        mkpath(output_dir)

        print_and_run_cmd(`fastqc --threads $(minimum((length(fastq_gzs), n_job))) --outdir $output_dir $fastq_gzs`)

        end_time = now()

        println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

    end

end
