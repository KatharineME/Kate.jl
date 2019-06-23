include("print_and_run_cmd.jl")


function find_mutation(
    germ_bam::String,
    soma_bam::String,
    dna_fa_bgz::String,
    chromosome_bed_gz::String,
    sequencing_scope::String,
    chrn_n_tsv::String,
    output_dir::String,
    n_job::Int,
    gb_memory::Int,
)

    println("Finding mutation ...")

    config_parameters::String = "--referencefasta $dna_fa_bgz --callregions $chromosome_bed_gz --$sequencing_scope"

    if ispath(germ_bam) && ispath(soma_bam)

        config_parameters = "$config_parameters --normalbam $germ_bam --tumorbam $soma_bam"

    elseif ispath(germ_bam)

        config_parameters = "$config_parameters --bam $germ_bam"

    else

        error("Arguments do not contain germ .bam.")

    end

    run_parameters::String = "--mode local --jobs $n_job --memGb $gb_memory --quiet"

    manta_dir::String = "$output_dir/manta"

    print_and_run_cmd(`bash -c "source activate py2.7 && configManta.py $config_parameters --outputContig --runDir $manta_dir && $manta_dir/runWorkflow.py $run_parameters"`)

    local configure_strelka::Cmd

    strelka_dir::String = "$output_dir/strelka"

    if ispath(germ_bam) && ispath(soma_bam)

        configure_strelka = `configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $manta_dir/results/variants/candidateSmallIndels.vcf.gz --runDir $strelka_dir`

    else

        configure_strelka = `configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_dir`

    end

    print_and_run_cmd(`bash -c "source activate py2.7 && $configure_strelka && $strelka_dir/runWorkflow.py $run_parameters"`)

    local concat_vcfs::String

    if ispath(germ_bam) && ispath(soma_bam)

        strelka_soma_sample_name_txt::String = "/tmp/strelka_soma_sample_name.txt"

        open(strelka_soma_sample_name_txt, "w") do io

            write(io, "Germ\nSoma")

        end

        indel_vcf_gz::String = "$strelka_dir/results/variants/somatic.indels.vcf.gz"

        tmp_indel_vcf_gz::String = "/tmp/indel.vcf.gz"

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $strelka_soma_sample_name_txt $indel_vcf_gz`,
            tmp_indel_vcf_gz,
        ))

        print_and_run_cmd(`mv --force $tmp_indel_vcf_gz $indel_vcf_gz`)

        print_and_run_cmd(`tabix --force $indel_vcf_gz`)

        snv_vcf_gz::String = "$strelka_dir/results/variants/somatic.snvs.vcf.gz"

        tmp_snv_vcf_gz::String = "/tmp/snv.vcf.gz"

        print_and_run_cmd(pipeline(
            `bcftools reheader --threads $n_job --samples $strelka_soma_sample_name_txt $snv_vcf_gz`,
            tmp_snv_vcf_gz,
        ))

        print_and_run_cmd(`mv --force $tmp_snv_vcf_gz $snv_vcf_gz`)

        print_and_run_cmd(`tabix --force $snv_vcf_gz`)

        print_and_run_cmd(`rm --force $strelka_soma_sample_name_txt`)

        concat_vcfs = "$manta_dir/results/variants/somaticSV.vcf.gz $indel_vcf_gz $snv_vcf_gz"

    else

        concat_vcfs = "$manta_dir/results/variants/diploidSV.vcf.gz $strelka_dir/results/variants/variants.vcf.gz"

    end

    print_and_run_cmd(pipeline(
        `bcftools concat --threads $n_job --allow-overlaps $concat_vcfs`,
        `bcftools annotate --threads $n_job --rename-chrs $chrn_n_tsv`,
        `bgzip --threads $n_job --stdout`,
        "/tmp/concat.vcf.gz",
    ))

    print_and_run_cmd(`tabix /tmp/concat.vcf.gz`)

    snpeff_dir::String = "$output_dir/snpeff"

    tmp_concat_vcf_gz::String = "/tmp/concat.vcf.gz"

    print_and_run_cmd(pipeline(
        `snpEff -Xms$(gb_memory//2)g -Xmx$(gb_memory)g GRCh38.86 -verbose -noLog -csvStats $snpeff_dir/stats.csv -htmlStats $snpeff_dir/stats.html $tmp_concat_vcf_gz`,
        `bgzip --threads $n_job --stdout`,
        "$snpeff_dir/variant.vcf.gz",
    ))

    print_and_run_cmd(`tabix $snpeff_dir/variant.vcf.gz`)

    print_and_run_cmd(`rm --force $tmp_concat_vcf_gz $tmp_concat_vcf_gz.tbi`)

end
