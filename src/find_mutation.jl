include("print_and_run_cmd.jl")


function find_mutation(
    germ_bam::String,
    soma_bam::String,
    dna_fa_bgz::String,
    chromosome_bed_gz::String,
    exome::Bool,
    chrn_n_tsv::String,
    output_dir::String,
    n_job::Int,
    gb_memory::Int,
)

    println("Finding mutation ...")

    config_parameters::String = "--referenceFasta $dna_fa_bgz --callRegions $chromosome_bed_gz"

    if exome

        config_parameters = "$config_parameters --exome"

    end

    if ispath(germ_bam) && ispath(soma_bam)

        config_parameters = "$config_parameters --normalbam $germ_bam --tumorbam $soma_bam"

    elseif ispath(germ_bam)

        config_parameters = "$config_parameters --bam $germ_bam"

    else

        error("Arguments do not contain germ .bam.")

    end

    run_parameters::String = "--mode local --jobs $n_job --memGb $gb_memory --quiet"

    partial_path::String = joinpath(
        "results",
        "variants",
    )

    manta_dir::String = joinpath(
        output_dir,
        "manta",
    )

    print_and_run_cmd(`bash -c "source activate py2.7 && configManta.py $config_parameters --outputContig --runDir $manta_dir && $(joinpath(manta_dir, "runWorkflow.py")) $run_parameters"`)

    strelka_dir::String = joinpath(
        output_dir,
        "strelka",
    )

    local configure_strelka::String

    if ispath(germ_bam) && ispath(soma_bam)

        candidatesmallindels_vcf_gz::String = joinpath(
            manta_dir,
            partial_path,
            "candidateSmallIndels.vcf.gz",
        )

        configure_strelka = "configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $candidatesmallindels_vcf_gz --runDir $strelka_dir"

    else

        configure_strelka = "configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_dir"

    end

    print_and_run_cmd(`bash -c "source activate py2.7 && $configure_strelka && $(joinpath(strelka_dir, "runWorkflow.py")) $run_parameters"`)

    local concat_vcfs::Tuple{Vararg{String}}

    if ispath(germ_bam) && ispath(soma_bam)

        sample_txt::String = joinpath(
            output_dir,
            "sample.txt",
        )

        open(sample_txt, "w") do io

            write(io, "Germ\nSoma")

        end

        indel_vcf_gz::String = joinpath(
            strelka_dir,
            partial_path,
            "somatic.indels.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $indel_vcf_gz`,
            "$indel_vcf_gz.tmp",
        ))

        print_and_run_cmd(`mv --force $indel_vcf_gz.tmp $indel_vcf_gz`)

        print_and_run_cmd(`tabix --force $indel_vcf_gz`)

        snv_vcf_gz::String = joinpath(
            strelka_dir,
            partial_path,
            "somatic.snvs.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $snv_vcf_gz`,
            "$snv_vcf_gz.tmp",
        ))

        print_and_run_cmd(`mv --force $snv_vcf_gz.tmp $snv_vcf_gz`)

        print_and_run_cmd(`tabix --force $snv_vcf_gz`)

        concat_vcfs = (
            joinpath(
                manta_dir,
                partial_path,
                "somaticSV.vcf.gz",
            ),
            indel_vcf_gz,
            snv_vcf_gz,
        )

    else

        concat_vcfs = (
            joinpath(
                manta_dir,
                partial_path,
                "diploidSV.vcf.gz",
            ),
            joinpath(
                strelka_dir,
                partial_path,
                "variants.vcf.gz",
            ),
        )

    end

    tmp_concat_vcf_gz = joinpath(
        output_dir,
        "concat.vcf.gz",
    )

    print_and_run_cmd(pipeline(
        `bcftools concat --threads $n_job --allow-overlaps $concat_vcfs`,
        `bcftools annotate --threads $n_job --rename-chrs $chrn_n_tsv`,
        `bgzip --threads $n_job --stdout`,
        tmp_concat_vcf_gz,
    ))

    print_and_run_cmd(`tabix $tmp_concat_vcf_gz`)

    snpeff_dir::String = joinpath(
        output_dir,
        "snpeff",
    )

    mkpath(snpeff_dir)

    variant_vcf_gz::String = joinpath(
        snpeff_dir,
        "variant.vcf.gz"
    )

    print_and_run_cmd(pipeline(
        `snpEff -Xmx$(gb_memory)g GRCh38.86 -verbose -noLog -csvStats $(joinpath(snpeff_dir, "stats.csv")) -htmlStats $(joinpath(snpeff_dir, "stats.html")) $tmp_concat_vcf_gz`,
        `bgzip --threads $n_job --stdout`,
        variant_vcf_gz,
    ))

    print_and_run_cmd(`tabix $variant_vcf_gz`)

    print_and_run_cmd(`rm --force $tmp_concat_vcf_gz $tmp_concat_vcf_gz.tbi`)

end
