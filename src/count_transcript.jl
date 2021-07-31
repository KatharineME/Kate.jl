using Dates

function count_transcript(
    fq1::String,
#     fq2::String,
    fa::String,
    pa::String,
    n_jo::Int,
)

    st = now()

    println("($st) Counting transcript ...")

    id::String = "$fa.kallisto_index"

    if !ispath(id)

        run_command(
            `kallisto index --index $id $fa`,
        )

    end

    mkpath(pa)

    run_command(
        `kallisto quant --single -l 51 -s 0.05 --threads $n_jo --index $id --output-dir $pa $fq1`,
    )
    
    # run_command(
    #    `kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`,
    #)


    en = now()

    println(
        "($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).",
    )

end

export count_transcript