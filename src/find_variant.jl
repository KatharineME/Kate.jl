using Dates

function find_variant(
    ge::Union{String,Nothing},
    so::Union{String,Nothing},
    ta::Bool,
    fa::String,
    chsi::String,
    chna::String,
    paou::String,
    n_jo::Int,
    me::Int,
    pasn::String,
)
    
    println("CALL VARIANTS")

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run_command(`samtools faidx $fa`)

    end
    
    if !ispath("$chsi.tbi")

        run_command(`tabix --force $chsi`)

    end

    # Set config parameters

    co::String = "--referenceFasta $fa --callRegions $chsi"

    if ta

        co = "$co --exome"

    end

    if ge != nothing && so != nothing

        co = "$co --normalBam $ge --tumorBam $so"

    elseif ge != nothing

        co = "$co --bam $ge"

    else

        error(
            "Specify germ and soma .bam for processing soma or germ .bam for processing germ.",
        )

    end

    # Set run parameters

    ru::String = "--mode local --jobs $n_jo --memGb $me --quiet"

    pava::String = joinpath("results", "variants")

    pama::String = joinpath(paou, "manta")

    run_command(
        `bash -c "source activate py2 && configManta.py $co --outputContig --runDir $pama && $(joinpath(pama, "runWorkflow.py")) $ru"`,
    )

    past::String = joinpath(paou, "strelka")

    # Configure strelka

    local st::String

    if ge != nothing && so != nothing

        st = "configureStrelkaSomaticWorkflow.py $co --indelCandidates $(joinpath(pama, pava, "candidateSmallIndels.vcf.gz")) --runDir $past"

    else

        st = "configureStrelkaGermlineWorkflow.py $co --runDir $past"

    end

    run_command(
        `bash -c "source activate py2 && $st && $(joinpath(past, "runWorkflow.py")) $ru"`,
    )

#     local vc_::Vector{Vararg{String}}

    if ge != nothing && so != nothing

        sa = joinpath(paou, "sample.txt")

        # TODO: get sample names (maybe from .bam) and use them instead of "Germ" and "Soma"

        open(io -> write(io, "Germ\nSoma"), sa; write = true)
        
        pain =
            joinpath(past, pava, "somatic.indels.vcf.gz")

        run_command(
            pipeline(
                `bcftools reheader --threads $n_jo --samples $sa $pain`,
                "$pain.tmp",
            ),
        )

        mv("$pain.tmp", pain; force = true)

        run_command(`tabix --force $pain`)

        pasv::String =
            joinpath(past, pava, "somatic.snvs.vcf.gz")

        run_command(
            pipeline(
                `bcftools reheader --threads $n_jo --samples $sa $pasv`,
                "pasv.tmp",
            ),
        )

        mv("$pasv.tmp", pasv; force = true)

        run_command(`tabix --force $pasv`)

        vc_ = [
            joinpath(pama, pava, "somaticSV.vcf.gz"),
            pain,
            pasv,
        ]

    else

        vc_ = [
            joinpath(pama, pava, "diploidSV.vcf.gz"),
            joinpath(past, pava, "variants.vcf.gz"),
           ]

    end
    
    paco::String = joinpath(paou, "concat.vcf.gz")

    run_command(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chna`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

    run_command(`tabix $paco`)

    sn::String = joinpath(paou, "snpeff")

    mkpath(sn)

    snvc::String = joinpath(sn, "snpeff.vcf.gz")

    run_command(
        pipeline(
            `java -Xmx$(me)g -jar $pasn GRCh38.99 -noLog -verbose -csvStats $(joinpath(sn, "stats.csv")) -htmlStats $(joinpath(sn, "stats.html")) $paco`,
            `bgzip --threads $n_jo --stdout`,
            snvc,
        ),
    )

    run_command(`tabix $snvc`)

    ps::String = joinpath(paou, "pass.vcf.gz")

    run_command(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $snvc`,
            `bgzip --threads $n_jo --stdout`,
            ps,
        ),
    )

    run_command(`tabix $ps`)
    
end

export find_variant