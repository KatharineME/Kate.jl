using Dates: now, CompoundPeriod

function concatenate_reads(fq_, sa::String, pa::String)

    st = now()

    fo_ = []

    re_ = []

    for fi in fq_

        if occursin("R1", fi)

            push!(fo_, fi)

        end
        
        if occursin("_1.fq", fi)
            
            push!(fo_, fi)
        
        end

    end

    for fi in fq_

        if occursin("R2", fi)

            push!(re_, fi)

        end
            
        if occursin("_2.fq", fi)
            
            push!(re_, fi)
    
        end

    end

    n_fo = length(fo_)

    n_re = length(re_)

    println("Number of forward (R1) read files = $n_fo\n")

    println("Number of reverse (R2) read files = $n_re\n")

    paca = joinpath(pa, string(sa, "_cat"))
        
    if ispath(paca)

        println(
            "Skipping concatenation because directory already exists:\n $paca\n",
        )
    
    elseif n_fo <= 1 && n_re <= 1

        println(
            "Nothing to concatenate because number of forward reads ($n_fo) and number of reverse reads ($n_re) are <= 1.",
                )
        
    else

        run(pipeline(`mkdir $paca`))        
            
        if n_fo > 1

            println("\nCombining R1 reads\n")

            run(
                pipeline(
                    `cat $fo_`,
                    stdout = joinpath(paca, string(sa, "_R1.fastq.gz")),
                ),
            )

        end

        if n_re > 1

            println("\nCombining R2 reads\n")

            run(
                pipeline(
                    `cat $re_`,
                    stdout = joinpath(
                        paca,
                        string(sa, "_R2.fastq.gz"),
                    ),
                ),
            )
        
        end

    end

    en = now()

    println("\nDone at: $en\n")

    println("Took $(canonicalize(Dates.CompoundPeriod(en - st))).\n")

end

export concatenate_reads