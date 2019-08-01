using Dates

include("print_and_run_cmd.jl")


function check_sequence(
    fastq_gz_file_paths::Tuple{Vararg{String}},
    output_directory_path::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Checking sequence ...")

    for fastq_gz_file_path in fastq_gz_file_paths

        html_file_path_suffix = "_fastqc.html"

        html_file_path = joinpath(
            output_directory_path,
            replace(
                replace(
                    splitdir(fastq_gz_file_path)[end],
                    ".fastq.gz" => html_file_path_suffix,
                ),
                ".fq.gz" => html_file_path_suffix,
            ),
        )

        if isfile(html_file_path)

            error("$html_file_path exists.")

        end

    end

    mkpath(output_directory_path)

    print_and_run_cmd(`fastqc --threads $n_job --outdir $output_directory_path $fastq_gz_file_paths`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
