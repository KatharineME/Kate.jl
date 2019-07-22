using Dates

include("print_and_run_cmd.jl")


function check_sequence(
    fastq_gzs::Tuple{Vararg{String}},
    output_dir::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Checking sequence ...")

    for fastq_gz in fastq_gzs

        suffix = "_fastqc.html"

        html = joinpath(
            output_dir,
            replace(
                replace(
                    split(
                        fastq_gz,
                        "/",
                    )[end],
                    ".fastq.gz"=>suffix,
                ),
                ".fq.gz"=>suffix,
            ),
        )

        println(html)

        if isfile(html)

            error("$html exists.")

        end

    end

    mkpath(output_dir)

    n_job_ = minimum((
        length(fastq_gzs),
        n_job,
    ))

    print_and_run_cmd(`fastqc --threads $n_job_ --outdir $output_dir $fastq_gzs`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
