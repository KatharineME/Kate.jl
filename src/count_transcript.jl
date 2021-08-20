using Dates

function count_transcript(
    fa::String,
    pa::String,
    n_jo::Int,
    fq1::String,
    fq2::String=nothing,
    fr::Int64=51,
    sd::Float64=0.05,  
)

    sa = now()

    println("Counting transcript ...\n")
    
    id::String = "$fa.kallisto_index"

    if !ispath(id)

        run_command(
            `kallisto index --index $id $fa`,
        )

    end

    mkpath(pa)
    
    if fq2 !== nothing

        run_command(
           `kallisto quant --threads $n_jo --index $id --output-dir $pa $fq1 $fq2`,
        )
        
    elseif fq2 === nothing
        
        run_command(
            `kallisto quant --single --fragment-length $fr --sd $sd --threads $n_jo --index $id --output-dir $pa $fq1`,
        )

    end
    
    en = now()

    println("Done at $en in $(canonicalize(Dates.CompoundPeriod(en - sa))).\n")

end

export count_transcript
