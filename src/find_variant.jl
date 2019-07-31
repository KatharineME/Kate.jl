using Dates

include("print_and_run_cmd.jl")


function find_variant(
    germ_bam::Union{
        String,
        Nothing,
    },
    soma_bam::Union{
        String,
        Nothing,
    },
    is_targeted::Bool,
    dna_fasta_bgz::String,
    chromosome_bed_gz::String,
    chrn_n_tsv::String,
    output_dir::String,
    n_job::Int,
    gb_memory::Int,
)

    start_time = now()

    println("($start_time) Finding variant ...")

    if !(isfile("$dna_fasta_bgz.fai") && ispath("$dna_fasta_bgz.gzi"))

        print_and_run_cmd(`samtools faidx $dna_fasta_bgz`)

    end

    if !ispath("$chromosome_bed_gz.tbi")

        print_and_run_cmd(`tabix --force $chromosome_bed_gz`)

    end

    config_parameters = `--referenceFasta $dna_fasta_bgz --callRegions $chromosome_bed_gz`

    if is_targeted

        config_parameters = `$config_parameters --exome`

    end

    # TODO: Check the best practice to check for nothing
    if germ_bam !== nothing && soma_bam !== nothing

        config_parameters = `$config_parameters --normalBam $germ_bam --tumorBam $soma_bam`

    elseif germ_bam !== nothing

        config_parameters = `$config_parameters --bam $germ_bam`

    else

        error("Specify germ_bam and soma_bam (for processing soma) or germ_bam (for processing germ).")

    end

    run_parameters = `--mode local --jobs $n_job --memGb $gb_memory --quiet`

    partial_path = joinpath(
        "results",
        "variants",
    )

    manta_dir = joinpath(
        output_dir,
        "manta",
    )

    print_and_run_cmd(`configManta.py $config_parameters --outputContig --runDir $manta_dir`)

    manta_runworkflow_py = joinpath(
        manta_dir,
        "runWorkflow.py",
    )

    print_and_run_cmd(`$manta_runworkflow_py $run_parameters`)

    strelka_dir = joinpath(
        output_dir,
        "strelka",
    )

    if germ_bam !== nothing && soma_bam !== nothing

        candidatesmallindels_vcf_gz = joinpath(
            manta_dir,
            partial_path,
            "candidateSmallIndels.vcf.gz",
        )

        print_and_run_cmd(`configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $candidatesmallindels_vcf_gz --runDir $strelka_dir`)

    else

        print_and_run_cmd(`configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_dir`)

    end

    strelka_runworkflow_py = joinpath(
        strelka_dir,
        "runWorkflow.py",
    )

    print_and_run_cmd(`$strelka_runworkflow_py $run_parameters`)

    if germ_bam !== nothing && soma_bam !== nothing

        sample_txt = joinpath(
            output_dir,
            "sample.txt",
        )

        # TODO: Use actual sample names instead of "Germ" and "Soma"
        open(
            io -> write(
                io,
                "Germ\nSoma",
            ),
            sample_txt;
            write = true,
        )

        somatic_indel_vcf_gz = joinpath(
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
           force = true,
        )

        print_and_run_cmd(`tabix --force $somatic_indel_vcf_gz`)

        somatic_snv_vcf_gz = joinpath(
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
           force = true,
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

    concat_vcf_gz = joinpath(
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

    snpeff_dir = joinpath(
        output_dir,
        "snpeff",
    )

    mkpath(snpeff_dir)

    snpeff_vcf_gz = joinpath(
        snpeff_dir,
        "snpeff.vcf.gz",
    )

    stats_csv = joinpath(
        snpeff_dir,
        "stats.csv",
    )

    stats_html = replace(
        stats_csv,
        ".csv" => ".html",
    )

    print_and_run_cmd(pipeline(
        `snpEff -Xmx$(gb_memory)g GRCh38.86 -noLog -verbose -csvStats $stats_csv -htmlStats $stats_html $concat_vcf_gz`,
        `bgzip --threads $n_job --stdout`,
        snpeff_vcf_gz,
    ))

    print_and_run_cmd(`tabix $snpeff_vcf_gz`)

    pass_vcf_gz = joinpath(
        output_dir,
        "pass.vcf.gz",
    )

    print_and_run_cmd(pipeline(
        `bcftools view --threads $n_job --include 'FILTER=="PASS"' $snpeff_vcf_gz`,
        `bgzip --threads $n_job --stdout`,
        pass_vcf_gz,
    ))

    print_and_run_cmd(`tabix $pass_vcf_gz`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
