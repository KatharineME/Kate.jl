using Dates

function count_transcript(
    fa::String,
    pa::String,
    n_jo::Int,
    fq1::String,
    fq2::String=nothing,
    st::String="pe",
    fr::Int64=51,
    sd::Float64=0.05,  
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
    
    if st == "pe"

        run_command(
           `kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`,
        )
        
    elseif st == "se"
        
        run_command(
            `kallisto quant --single --fragment-length $fr --sd $sd --threads $n_jo --index $id --output-dir $pa $fq1`,
        )

    end
    
    en = now()

    println("Done at $en in $(canonicalize(Dates.CompoundPeriod(en - st))).\n")

end

export count_transcript