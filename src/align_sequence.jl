using Dates

include("run_command.jl")


function align_sequence(
    fq1::String,
    fq2::String,
    sa::String,
    fa::String,
    pa::String,
    n_jo::Int,
    me::Int,
)

    st = now()

    println("($st) Aligning sequence ...")

    id::String = "$fa.mmi"

    if !ispath(id)

        # Make index
        #
        run_command(`minimap2 -t $n_jo -d $id $fa`)

    end
    
    paal = splitdir(pa)[1]
    
    if isdir(paal)
         
        continue
        
    else
        
        mkpath(paal)

        run_command(
            pipeline(
                `minimap2 -x sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" -a $id $fq1 $fq2`,
                # `samtools sort --threads $n_jo -m $(me)G -n`,
                `samtools sort --threads $n_jo -n`,
                `samtools fixmate --threads $n_jo -m - -`,
                # `samtools sort --threads $n_jo -m $(me)G`,
                `samtools sort --threads $n_jo`,
                "$pa.tmp",
            ),
        )

        run_command(`samtools markdup --threads $n_jo -s $pa.tmp $pa`)

        rm("$pa.tmp")

        run_command(`samtools index -@ $n_jo $pa`)

        run_command(pipeline(`samtools flagstat --threads $n_jo $pa`, "$pa.flagstat"))

    end
        
    en = now()

    println(
        "($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).",
    )

end
