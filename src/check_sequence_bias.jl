using Dates

include("run_command.jl")


function check_sequence_bias(pa::String)

    st = now()

    run_command(
        `multiqc --outdir $pa $pa`,
    )

    en = now()

    println("\nDone at: $en\n")

    println("Took $(canonicalize(Dates.CompoundPeriod(en - st))).\n")

end
