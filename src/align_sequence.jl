using Dates: now, CompoundPeriod

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

    println("ALIGN")
    
    paal = splitdir(pa)[1]
    
    id::String = "$fa.mmi"
    
    if isdir(paal)
         
        println("Skipping alignment because directory already exists: \n$paal")
        
    else
        
        println("($st) Aligning sequence ...")
        
        if !ispath(id)

            run_command(`minimap2 -t $n_jo -d $id $fa`)

        end
        
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
        "($en) Done in $(canonicalize(Dates.CompoundPeriod(en - st))).\n",
    )

end

export align_sequence