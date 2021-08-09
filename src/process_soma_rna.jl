function process_soma_rna(
    so1::String,
    so2::String,
    pa::String,
    fa::String,  # fa is cdna_fasta.gz
    n_jo::Int,
)

    for fi in (so1, so2, fa)
        if !isfile(fi)

            error("$fi doesn't exist.")

        end

    end

    patr = joinpath(pa, "trim_sequence", "soma")

    trim_sequence(
        so1,
        so2,
        patr,
        n_jo,
    )

    so1tr = "$patr-trimmed-pair1.fastq.gz"

    so2tr = "$patr-trimmed-pair2.fastq.gz"

    check_sequence(
        (so1tr, so2tr),
        joinpath(pa, "check_sequence"),
        n_jo,
    )

    paco = joinpath(pa, "count_transcript")

    count_transcript(
        so1tr,
        so2tr,
        fa,
        paco,
        n_jo,
    )

end

export process_soma_rna
