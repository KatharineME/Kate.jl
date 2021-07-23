using Dates

include("run_command.jl")


function check_sequence(fq_::Array, pa::String, n_jo::Int)

    st = now()
    
    if ispath(pa)

        println(
            "Skipping check sequence because directory already exists:\n $pa\n",
        )

    else

        println("($st) Checking sequence ...")

        mkpath(pa)

        run_command(
            `fastqc --threads $(minimum((length(fq_), n_jo))) --outdir $pa $fq_`,
        )

        en = now()

        println(
            "($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).",
        )

    end

end
