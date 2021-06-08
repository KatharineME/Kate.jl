using Dates

include("run_command.jl")


function count_transcript(
    fq1::String,
    fq2::String,
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
        `kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`,
    )

    en = now()

    println(
        "($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).",
    )

end
