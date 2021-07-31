using Dates

function find_read(pa::String)

    st = now()

    n_fq = 0

    n_gz = 0

    fi_ = []

    println("Walking sample directory...\n")

    for (ro, di, th_) in walkdir("$pa")

        println("$ro\n")

        for th in th_

            if occursin(".fastq", th) && !(occursin(".md5", th))

                n_fq += 1

            end

            if occursin("fastq.gz", th) && !(occursin(".md5", th))

                n_gz += 1

                push!(fi_, joinpath(ro, th))

            end

            if occursin("fq.gz", th) && !(occursin(".md5", th))

                n_gz += 1

                push!(fi_, joinpath(ro, th))

            end

        end

    end

    println("\nNumber of fastq files found in directories walked: $n_fq\n")

    println(
        "Number of fastq.gz or fq.gz files found in directories walked: $n_gz\n",
    )

    en = now()

    println("\nDone at: $en\n")

    println("Took $(canonicalize(Dates.CompoundPeriod(en - st))).\n")

    return fi_

end

export find_read