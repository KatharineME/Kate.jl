using Dates

include("run_command.jl")


function check_sequence_bias(sa::String, pa::String)

    st = now()

    pach = joinpath(pa, string("check_sequence_", sa))

    run_command(
        `multiqc --outdir $pach $pach`,
    )

    en = now()

    println("\nDone at: $en\n")

    println("Took $(canonicalize(Dates.CompoundPeriod(en - st))).\n")

end
