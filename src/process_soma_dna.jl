include("trim_sequence.jl")

include("check_sequence.jl")

include("align_sequence.jl")

include("find_variant.jl")


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

    germ_trim_sequence_prefix = joinpath(output_dir, "trim_sequence", "germ")

    trim_sequence(
        germ_dna_1_fastq_gz,
        germ_dna_2_fastq_gz,
        germ_trim_sequence_prefix,
        n_job,
    )

    germ_trim_1_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair1.fastq.gz"

    germ_trim_2_fastq_gz = "$germ_trim_sequence_prefix-trimmed-pair2.fastq.gz"

    soma_trim_sequence_prefix = joinpath(output_dir, "trim_sequence", "soma")

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
        joinpath(output_dir, "check_sequence"),
        n_job,
    )

    germ_bam = joinpath(output_dir, "align_sequence", "germ.bam")

    align_sequence(
        germ_trim_1_fastq_gz,
        germ_trim_2_fastq_gz,
        "Germ",
        dna_fasta_gz,
        germ_bam,
        n_job,
        job_gb_memory,
    )

    soma_bam = joinpath(output_dir, "align_sequence", "soma.bam")

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

    find_variant_dir = joinpath(output_dir, "find_variant")

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
