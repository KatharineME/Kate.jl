using Dates: now, CompoundPeriod

function trim_sequence(fq1::String, fq2::String, pa::String, n_jo::Int)

    st = now()
    
    println("TRIM")

    pa1 = joinpath(string(pa), "trimmed-pair1.fastq.gz")

    pa2 = joinpath(string(pa), "trimmed-pair2.fastq.gz")

    if isfile(pa1) && isfile(pa2)

        println("Skipping trimming because trimmed files already exist:\n $pa1\n $pa2\n")

    else


        println("($st) Trimming sequence ...")

        mkpath(pa)
        
        println("Made path for trimmed files: $pa")

        run_command(
            `skewer --threads $n_jo -x AGATCGGAAGAGC --mode pe --mean-quality 10 --end-quality 15 --compress --output $pa --quiet $fq1 $fq2`,
        )

        en = now()

        println("($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).\n\n")

    end

end

export trim_sequence
