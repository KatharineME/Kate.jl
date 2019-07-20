module Kraft

using Dates


function print_and_run_cmd(
    cmd::Base.AbstractCmd,
)

    println(cmd)

    run(cmd)

end


function error_if_same(
    _1::String,
    _2::String,
)

    if _1 == _2

        error("$_1 and $_2 are the same.")

    end

end


function check_sequence(
    fastq_gzs::Tuple{Vararg{String}},
    output_dir::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Checking sequence ...")

    for fastq_gz in fastq_gzs

        suffix = "_fastqc.html"

        html = joinpath(
            output_dir,
            replace(
                replace(
                    split(
                        fastq_gz,
                        "/",
                    )[end],
                    ".fastq.gz"=>suffix,
                ),
                ".fq.gz"=>suffix,
            ),
        )

        println(html)

        if isfile(html)

            error("$html exist.")

        end

    end

    mkpath(output_dir)

    print_and_run_cmd(`fastqc --threads $(minimum((
        length(fastq_gzs),
        n_job,
    ))) --outdir $output_dir $fastq_gzs`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end


function trim_sequence(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    output_prefix::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Trimming sequence ...")

    if isfile("$output_prefix-trimmed-pair1.fastq.gz") || isfile("$output_prefix-trimmed-pair2.fastq.gz")

        error("$output_prefix-trimmed-pair(1|2).fastq.gz exists.")

    end
    
    output_dir = splitdir(output_prefix)[1]

    mkpath(output_dir)

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_prefix --quiet $_1_fastq_gz $_2_fastq_gz`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end


function align_sequence(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    sample_name::String,
    dna_fasta_gz::String,
    bam::String,
    n_job::Int,
    job_gb_memory::Int,
)

    start_time = now()

    println("($start_time) Aligning sequence ...")

    dna_fasta_gz_mmi = "$dna_fasta_gz.mmi"

    if !ispath(dna_fasta_gz_mmi)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fasta_gz_mmi $dna_fasta_gz`)

    end

    if isfile(bam)

        error("$bam exists.")
    
    end

    output_dir = splitdir(bam)[1]

    mkpath(output_dir)

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -K $(job_gb_memory)G -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fasta_gz_mmi $_1_fastq_gz $_2_fastq_gz`,
        `samtools sort --threads $n_job -m $(job_gb_memory)G -n`,
        `samtools fixmate --threads $n_job -m - -`,
        `samtools sort --threads $n_job -m $(job_gb_memory)G`,
        "$bam.tmp",
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $bam.tmp $bam`)

    rm("$bam.tmp")

    print_and_run_cmd(`samtools index -@ $n_job $bam`)

    print_and_run_cmd(pipeline(
        `samtools flagstat --threads $n_job $bam`,
        "$bam.flagstat",
    ))

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end


function find_variant(
    germ_bam::Union{String, Nothing},
    soma_bam::Union{String, Nothing},
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

    if germ_bam != nothing && soma_bam != nothing

        config_parameters = `$config_parameters --normalBam $germ_bam --tumorBam $soma_bam`

    elseif germ_bam != nothing

        config_parameters = `$config_parameters --bam $germ_bam`

    else

        error("Specify germ and soma .bam for processing soma or germ .bam for processing germ.")

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

    print_and_run_cmd(`$(joinpath(
        manta_dir,
        "runWorkflow.py",
    )) $run_parameters`)

    strelka_dir = joinpath(
        output_dir,
        "strelka",
    )

    if germ_bam != nothing && soma_bam != nothing

        print_and_run_cmd(`configureStrelkaSomaticWorkflow.py $config_parameters --indelCandidates $(joinpath(
            manta_dir,
            partial_path,
            "candidateSmallIndels.vcf.gz",
        )) --runDir $strelka_dir`)

    else

        print_and_run_cmd(`configureStrelkaGermlineWorkflow.py $config_parameters --runDir $strelka_dir`)

    end

    print_and_run_cmd(`$(joinpath(
        strelka_dir,
        "runWorkflow.py",
    )) $run_parameters`)

    if germ_bam != nothing && soma_bam != nothing

        sample_txt = joinpath(
            output_dir,
            "sample.txt",
        )

        # TODO: Use actual sample names instead of "Germ" and "Soma"
        open(
            io->write(
                io,
                "Germ\nSoma",
            ),
            sample_txt;
            write=true,
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
            force=true,
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

    print_and_run_cmd(pipeline(
        `snpEff -Xmx$(gb_memory)g GRCh38.86 -noLog -verbose -csvStats $(joinpath(
            snpeff_dir,
            "stats.csv",
        )) -htmlStats $(joinpath(
            snpeff_dir,
            "stats.html",
        )) $concat_vcf_gz`,
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

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end


function count_transcript(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    cdna_fasta_gz::String,
    output_dir::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Counting transcript ...")

    cdna_fasta_gz_kallisto_index = "$cdna_fasta_gz.kallisto_index"

    if !ispath(cdna_fasta_gz_kallisto_index)

        print_and_run_cmd(`kallisto index --index $cdna_fasta_gz_kallisto_index $cdna_fasta_gz`)

    end

    if isdir(output_dir)

        error("$output_dir exists.")
    
    else

        mkpath(output_dir)
    
    end

    print_and_run_cmd(`kallisto quant --threads $n_job --index $cdna_fasta_gz_kallisto_index --output-dir $output_dir $_1_fastq_gz $_2_fastq_gz`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end


function process_germ_dna(
    germ_dna_1_fastq_gz::String,
    germ_dna_2_fastq_gz::String,
    dna_is_targeted::Bool,
    output_dir::String,
    dna_fasta_gz::String,
    chromosome_bed_gz::String,
    chrn_n_tsv::String,
    n_job::Int,
    gb_memory::Int,
    job_gb_memory::Int,
)

    for file_path in (
        germ_dna_1_fastq_gz,
        germ_dna_2_fastq_gz,
        dna_fasta_gz,
        chromosome_bed_gz,
        chrn_n_tsv,
    )

        if !isfile(file_path)

            error("$file_path does not exist.")

        end

    end

    if isdir(output_dir)

        # error("$output_dir exists.")
    
    else

        mkpath(output_dir)
    
    end

    germ_trim_sequence_prefix = joinpath(
        output_dir,
        "trim_sequence",
        "germ",
    )

    # trim_sequence(
    #     germ_dna_1_fastq_gz,
    #     germ_dna_2_fastq_gz,
    #     germ_trim_sequence_prefix,
    #     n_job,
    # )

    germ_trim_1_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    germ_trim_2_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    # check_sequence(
    #     (
    #         germ_trim_1_fastq_gz,
    #         germ_trim_2_fastq_gz,
    #     ),
    #     joinpath(
    #         output_dir,
    #         "check_sequence",
    #     ),
    #     n_job,
    # )

    germ_bam = joinpath(
        output_dir,
        "align_sequence",
        "germ.bam",
    )

    # align_sequence(
    #     germ_trim_1_fastq_gz,
    #     germ_trim_2_fastq_gz,
    #     "Germ",
    #     dna_fasta_gz,
    #     germ_bam,
    #     n_job,
    #     job_gb_memory,
    # )

    dna_fasta_bgz = "$(splitext(dna_fasta_gz)[1]).bgz"

    if !isfile(dna_fasta_bgz)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fasta_gz --stdout`,
            `bgzip --threads $n_job --stdout`,
            dna_fasta_bgz,
        ))

    end

    find_variant_dir = joinpath(
        output_dir,
        "find_variant",
    )

    find_variant(
        germ_bam,
        nothing,
        dna_is_targeted,
        dna_fasta_bgz,
        chromosome_bed_gz,
        chrn_n_tsv,
        find_variant_dir,
        n_job,
        gb_memory,
    )

end


function process_soma_dna(
    germ_dna_1_fastq_gz::String,
    germ_dna_2_fastq_gz::String,
    soma_dna_1_fastq_gz::String,
    soma_dna_2_fastq_gz::String,
    dna_is_targeted::Bool,
    output_dir::String,
    dna_fasta_gz::String,
    chromosome_bed_gz::String,
    chrn_n_tsv::String,
    n_job::Int,
    gb_memory::Int,
    job_gb_memory::Int,
)

    for file_path in (
        germ_dna_1_fastq_gz,
        germ_dna_2_fastq_gz,
        soma_dna_1_fastq_gz,
        soma_dna_2_fastq_gz,
        dna_fasta_gz,
        chromosome_bed_gz,
        chrn_n_tsv,
    )

        if !isfile(file_path)

            error("$file_path does not exist.")

        end

    end

    if isdir(output_dir)

        error("$output_dir exists.")
    
    else

        mkpath(output_dir)
    
    end

    germ_trim_sequence_prefix = joinpath(
        output_dir,
        "trim_sequence",
        "germ",
    )

    trim_sequence(
        germ_dna_1_fastq_gz,
        germ_dna_2_fastq_gz,
        germ_trim_sequence_prefix,
        n_job,
    )

    germ_trim_1_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    germ_trim_2_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    soma_trim_sequence_prefix = joinpath(
        output_dir,
        "trim_sequence",
        "soma",
    )

    trim_sequence(
        soma_dna_1_fastq_gz,
        soma_dna_2_fastq_gz,
        soma_trim_sequence_prefix,
        n_job,
    )

    soma_trim_1_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    soma_trim_2_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    check_sequence(
        (
            germ_trim_1_fastq_gz,
            germ_trim_2_fastq_gz,
            soma_trim_1_fastq_gz,
            soma_trim_2_fastq_gz,
        ),
        joinpath(
            output_dir,
            "check_sequence",
        ),
        n_job,
    )

    germ_bam = joinpath(
        output_dir,
        "align_sequence",
        "germ.bam",
    )

    align_sequence(
        germ_trim_1_fastq_gz,
        germ_trim_2_fastq_gz,
        "Germ",
        dna_fasta_gz,
        germ_bam,
        n_job,
        job_gb_memory,
    )

    soma_bam = joinpath(
        output_dir,
        "align_sequence",
        "soma.bam",
    )

    align_sequence(
        soma_trim_1_fastq_gz,
        soma_trim_2_fastq_gz,
        "Soma",
        dna_fasta_gz,
        soma_bam,
        n_job,
        job_gb_memory,
    )

    dna_fasta_bgz = "$(splitext(dna_fasta_gz)[1]).bgz"

    if !isfile(dna_fasta_bgz)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fasta_gz --stdout`,
            `bgzip --threads $n_job --stdout`,
            dna_fasta_bgz,
        ))

    end

    find_variant_dir = joinpath(
        output_dir,
        "find_variant",
    )

    find_variant(
        germ_bam,
        soma_bam,
        dna_is_targeted,
        dna_fasta_bgz,
        chromosome_bed_gz,
        chrn_n_tsv,
        find_variant_dir,
        n_job,
        gb_memory,
    )

end


function process_soma_rna(
    soma_rna_1_fastq_gz::String,
    soma_rna_2_fastq_gz::String,
    output_dir::String,
    cdna_fasta_gz::String,
    n_job::Int,
)

    for file_path in (
        soma_rna_1_fastq_gz,
        soma_rna_2_fastq_gz,
        cdna_fasta_gz,
    )

        if !isfile(file_path)

            error("$file_path does not exist.")

        end

    end

    if isdir(output_dir)

        error("$output_dir exists.")
    
    else

        mkpath(output_dir)
    
    end

    soma_trim_sequence_prefix = joinpath(
        output_dir,
        "trim_sequence",
        "soma",
    )

    trim_sequence(
        soma_rna_1_fastq_gz,
        soma_rna_2_fastq_gz,
        soma_trim_sequence_prefix,
        n_job,
    )

    soma_trim_1_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    soma_trim_2_fastq_gz = "$soma_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    check_sequence(
        (
            soma_trim_1_fastq_gz,
            soma_trim_2_fastq_gz,
        ),
        joinpath(
            output_dir,
            "check_sequence",
        ),
        n_job,
    )

    count_transcript_dir = joinpath(
        output_dir,
        "count_transcript",
    )

    count_transcript(
        soma_trim_1_fastq_gz,
        soma_trim_2_fastq_gz,
        cdna_fasta_gz,
        count_transcript_dir,
        n_job,
    )

end


end
