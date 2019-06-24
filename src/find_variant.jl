using Dates

include("print_and_run_cmd.jl")


function find_variant(
    germ_bam::Union{String, Nothing},
    soma_bam::Union{String, Nothing},
    dna_fa_bgz::String,
    chromosome_bed_gz::String,
    is_targeted::Bool,
    chrn_n_tsv::String,
    output_dir::String,
    n_job::Int,
    gb_memory::Int,
)

    start_time = now()

    println("($start_time) Finding variant ...")

    if !ispath("$chromosome_bed_gz.tbi")

        print_and_run_cmd(`tabix --force $chromosome_bed_gz`)

    end

    config_parameters::String = "--referenceFasta $dna_fa_bgz --callRegions $chromosome_bed_gz"

    if is_targeted

        config_parameters = "$config_parameters --exome"

    end

    if germ_bam != nothing && soma_bam != nothing

        config_parameters = "$config_parameters --normalBam $germ_bam --tumorBam $soma_bam"

    elseif germ_bam != nothing

        config_parameters = "$config_parameters --bam $germ_bam"

    else

        error("Specify germ and soma .bam for processing soma or germ .bam for processing germ.")

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

    if germ_bam != nothing && soma_bam != nothing

        configure_strelka = "configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $(joinpath(manta_dir, partial_path, "candidateSmallIndels.vcf.gz")) --runDir $strelka_dir"

    else

        configure_strelka = "configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_dir"

    end

    print_and_run_cmd(`bash -c "source activate py2.7 && $configure_strelka && $(joinpath(strelka_dir, "runWorkflow.py")) $run_parameters"`)

    local concat_vcfs::Tuple{Vararg{String}}

    if germ_bam != nothing && soma_bam != nothing

        sample_txt::String = joinpath(
            output_dir,
            "sample.txt",
        )

        #TODO: get sample names (maybe from .bam) and use them instead of "Germ" and "Soma"
        open(
            io->write(io, "Germ\nSoma"),
            sample_txt;
            write=true,
        )

        somatic_indel_vcf_gz::String = joinpath(
            strelka_dir,
            partial_path,
            "somatic.indels.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $somatic_indel_vcf_gz`,
            "$somatic_indel_vcf_gz.tmp",
        ))

        mv(
            "$somatic_indel_vcf_gz.tmp",
            somatic_indel_vcf_gz;
            force=true,
        )

        print_and_run_cmd(`tabix --force $somatic_indel_vcf_gz`)

        somatic_snv_vcf_gz::String = joinpath(
            strelka_dir,
            partial_path,
            "somatic.snvs.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $somatic_snv_vcf_gz`,
            "$somatic_snv_vcf_gz.tmp",
        ))

        mv(
            "$somatic_snv_vcf_gz.tmp",
            somatic_snv_vcf_gz;
            force=true,
        )

        print_and_run_cmd(`tabix --force $somatic_snv_vcf_gz`)

        concat_vcfs = (
            joinpath(
                manta_dir,
                partial_path,
                "somaticSV.vcf.gz",
            ),
            somatic_indel_vcf_gz,
            somatic_snv_vcf_gz,
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

    concat_vcf_gz::String = joinpath(
        output_dir,
        "concat.vcf.gz",
    )

    print_and_run_cmd(pipeline(
        `bcftools concat --threads $n_job --allow-overlaps $concat_vcfs`,
        `bcftools annotate --threads $n_job --rename-chrs $chrn_n_tsv`,
        `bgzip --threads $n_job --stdout`,
        concat_vcf_gz,
    ))

    print_and_run_cmd(`tabix $concat_vcf_gz`)

    snpeff_dir::String = joinpath(
        output_dir,
        "snpeff",
    )

    mkpath(snpeff_dir)

    snpeff_vcf_gz::String = joinpath(
        snpeff_dir,
        "snpeff.vcf.gz",
    )

    print_and_run_cmd(pipeline(
        `snpEff -Xmx$(gb_memory)g GRCh38.86 -noLog -verbose -csvStats $(joinpath(snpeff_dir, "stats.csv")) -htmlStats $(joinpath(snpeff_dir, "stats.html")) $concat_vcf_gz`,
        `bgzip --threads $n_job --stdout`,
        snpeff_vcf_gz,
    ))

    print_and_run_cmd(`tabix $snpeff_vcf_gz`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end
