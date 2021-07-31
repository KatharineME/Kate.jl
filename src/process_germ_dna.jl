function process_germ_dna(
    fq1::String,
    fq2::String,
    ta::Bool,
    paou::String,
    fa::String,
    chsi::String,
    chna::String,
    n_jo::Int,
    meto::Int,
    mejo::Int,
    pasn::String,
)

    for file_path::String in (fq1, fq2, fa, chsi, chna, pasn)
        if !isfile(file_path)

            error("$file_path doesn't exist.")

        end

    end

    patr::String = joinpath(paou, "trim_sequence/")
    
    trim_sequence(fq1, fq2, patr, n_jo)

    fq1tr::String = "$patr/trimmed-pair1.fastq.gz"

    fq2tr::String = "$patr/trimmed-pair2.fastq.gz"

    check_sequence([fq1tr, fq2tr], joinpath(paou, "check_sequence_trimmed"), n_jo)

    paal::String = joinpath(paou, "align_sequence", "germ.bam")

    align_sequence(fq1tr, fq2tr, "Germ", fa, paal, n_jo, mejo)

    fagz::String = "$(splitext(fa)[1]).bgz"

    if !isfile(fagz)

        run_command(
            pipeline(
                `gzip --decompress $fa --stdout`,
                `bgzip --threads $n_jo --stdout`,
                fagz,
            ),
        )

    end

    pava::String = joinpath(paou, "find_variant")
    
    find_variant(paal, nothing, ta, fagz, chsi, chna, pava, n_jo, meto, pasn)

end

export process_germ_dna