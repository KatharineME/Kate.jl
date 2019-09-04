using Dates

include("print_and_run_cmd.jl")


function find_variant(
    germ_bam_file_path::Union{Nothing,String,},
    soma_bam_file_path::Union{Nothing,String,},
    is_targeted::Bool,
    dna_fasta_bgz_file_path::String,
    chromosome_bed_gz_file_path::String,
    chrn_n_tsv_file_path::String,
    output_directory_path::String,
    n_job::Int,
    n_gb_memory::Int,
)

    start_time = now()

    println("($start_time) Finding variant...")

    if !(isfile("$dna_fasta_bgz_file_path.fai") && isfile("$dna_fasta_bgz_file_path.gzi"))

        print_and_run_cmd(`samtools faidx $dna_fasta_bgz_file_path`)

    end

    if !isfile("$chromosome_bed_gz_file_path.tbi")

        print_and_run_cmd(`tabix --force $chromosome_bed_gz_file_path`)

    end

    config_parameters = `--referenceFasta $dna_fasta_bgz_file_path --callRegions $chromosome_bed_gz_file_path`

    if is_targeted

        config_parameters = `$config_parameters --exome`

    end

    if germ_bam_file_path !== nothing && soma_bam_file_path !== nothing

        config_parameters = `$config_parameters --normalBam $germ_bam_file_path --tumorBam $soma_bam_file_path`

    elseif germ_bam_file_path !== nothing

        config_parameters = `$config_parameters --bam $germ_bam_file_path`

    else

        error("germ_bam_file_path and soma_bam_file_path (for processing soma) or germ_bam_file_path (for processing germ) is missing.")

    end

    run_parameters = `--mode local --jobs $n_job --memGb $n_gb_memory --quiet`

    partial_path = joinpath("results", "variants",)

    manta_directory_path = joinpath(output_directory_path, "manta",)

    print_and_run_cmd(`configManta.py $config_parameters --outputContig --runDir $manta_directory_path`)

    manta_runworkflow_py = joinpath(manta_directory_path, "runWorkflow.py",)

    print_and_run_cmd(`$manta_runworkflow_py $run_parameters`)

    strelka_directory_path = joinpath(output_directory_path, "strelka",)

    if germ_bam_file_path !== nothing && soma_bam_file_path !== nothing

        candidatesmallindels_vcf_gz_file_path = joinpath(
            manta_directory_path,
            partial_path,
            "candidateSmallIndels.vcf.gz",
        )

        print_and_run_cmd(`configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $candidatesmallindels_vcf_gz_file_path --runDir $strelka_directory_path`)

    else

        print_and_run_cmd(`configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_directory_path`)

    end

    strelka_runworkflow_py = joinpath(strelka_directory_path, "runWorkflow.py",)

    print_and_run_cmd(`$strelka_runworkflow_py $run_parameters`)

    if germ_bam_file_path !== nothing && soma_bam_file_path !== nothing

        sample_txt = joinpath(output_directory_path, "sample.txt",)

        open(io -> write(io, "Germ\nSoma",), sample_txt; write = true,)

        somatic_indel_vcf_gz_file_path = joinpath(
            strelka_directory_path,
            partial_path,
            "somatic.indels.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $somatic_indel_vcf_gz_file_path`,
            "$somatic_indel_vcf_gz_file_path.tmp",
        ))

        mv(
           "$somatic_indel_vcf_gz_file_path.tmp",
           somatic_indel_vcf_gz_file_path;
           force = true,
        )

        print_and_run_cmd(`tabix --force $somatic_indel_vcf_gz_file_path`)

        somatic_snv_vcf_gz_file_path = joinpath(
            strelka_directory_path,
            partial_path,
            "somatic.snvs.vcf.gz",
        )

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $sample_txt $somatic_snv_vcf_gz_file_path`,
            "$somatic_snv_vcf_gz_file_path.tmp",
        ))

        mv("$somatic_snv_vcf_gz_file_path.tmp", somatic_snv_vcf_gz_file_path; force = true,)

        print_and_run_cmd(`tabix --force $somatic_snv_vcf_gz_file_path`)

        concat_vcf_file_paths = (
            joinpath(manta_directory_path, partial_path, "somaticSV.vcf.gz",),
            somatic_indel_vcf_gz_file_path,
            somatic_snv_vcf_gz_file_path,
        )

    else

        concat_vcf_file_paths = (
            joinpath(manta_directory_path, partial_path, "diploidSV.vcf.gz",),
            joinpath(strelka_directory_path, partial_path, "variants.vcf.gz",),
        )

    end

    concat_vcf_gz_file_path = joinpath(output_directory_path, "concat.vcf.gz",)

    print_and_run_cmd(pipeline(
        `bcftools concat --threads $n_job --allow-overlaps $concat_vcf_file_paths`,
        `bcftools annotate --threads $n_job --rename-chrs $chrn_n_tsv_file_path`,
        `bgzip --threads $n_job --stdout`,
        concat_vcf_gz_file_path,
    ))

    print_and_run_cmd(`tabix $concat_vcf_gz_file_path`)

    snpeff_directory_path = joinpath(output_directory_path, "snpeff",)

    mkpath(snpeff_directory_path)

    snpeff_vcf_gz_file_path = joinpath(snpeff_directory_path, "snpeff.vcf.gz",)

    stats_csv = joinpath(snpeff_directory_path, "stats.csv",)

    stats_html = replace(stats_csv, ".csv" => ".html",)

    print_and_run_cmd(pipeline(
        `snpEff -Xmx$(n_gb_memory)g GRCh38.86 -noLog -verbose -csvStats $stats_csv -htmlStats $stats_html $concat_vcf_gz_file_path`,
        `bgzip --threads $n_job --stdout`,
        snpeff_vcf_gz_file_path,
    ))

    print_and_run_cmd(`tabix $snpeff_vcf_gz_file_path`)

    pass_vcf_gz_file_path = joinpath(output_directory_path, "pass.vcf.gz",)

    print_and_run_cmd(pipeline(
        `bcftools view --threads $n_job --include 'FILTER=="PASS"' $snpeff_vcf_gz_file_path`,
        `bgzip --threads $n_job --stdout`,
        pass_vcf_gz_file_path,
    ))

    print_and_run_cmd(`tabix $pass_vcf_gz_file_path`)

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

    return nothing

end
