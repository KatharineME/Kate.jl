include("print_and_run_cmd.jl")


function make_path_dict(
    data_dir::String,
    remake::Bool,
    n_job::Int,
)

    println("Making path dict ...")

    dna_fa_gz::String = joinpath(
        data_dir,
        "grch",
        "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    )

    chromosome_bed_gz::String = joinpath(
        data_dir,
        "grch",
        "chromosome.bed.gz",
    )

    chrn_n_tsv::String = joinpath(
        data_dir,
        "grch",
        "chrn_n.tsv",
    )

    cdna_fa_gz::String = joinpath(
        data_dir,
        "grch",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
    )

    grch_enst_gene_name_tsv::String = joinpath(
        data_dir,
        "grch",
        "enst_gene_name.tsv",
    )

    virus_cdna_fa_gz::String = joinpath(
        data_dir,
        "virus",
        "sequences.fasta",
    )

    virus_sequences_csv::String = replace(
        virus_cdna_fa_gz,
        ".tsv"=>".csv",
    )

    for path::String in (
        dna_fa_gz,
        chromosome_bed_gz,
        chrn_n_tsv,
        cdna_fa_gz,
        grch_enst_gene_name_tsv,
        virus_cdna_fa_gz,
        virus_sequences_csv,
    )

        if ispath(path)

            println(path)

        else

            error("$path doesn't exist.")

        end

    end

    dna_fa_bgz::String = replace(
        dna_fa_gz,
        ".gz"=>".bgz",
    )

    if remake || !ispath(dna_fa_bgz)

        print_and_run_cmd(pipeline(
            `gzip --decompress $dna_fa_gz --stdout`,
            `bgzip --threads $n_job --stdout`,
            dna_fa_bgz,
        ))

    end

    dna_fa_gz_mmi = "$dna_fa_gz.mmi"

    if remake || !ispath(dna_fa_gz_mmi)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fa_gz_mmi $dna_fa_gz`)

    end

    if remake || !(ispath("$dna_fa_bgz.fai") && ispath("$dna_fa_bgz.gzi"))

        print_and_run_cmd(`samtools faidx $dna_fa_bgz`)

    end

    if remake || !ispath("$chromosome_bed_gz.tbi")

        print_and_run_cmd(`tabix --force $chromosome_bed_gz`)

    end

    cdna_fa_gz_kallisto_index = "$cdna_fa_gz.kallisto_index"

    if remake || !ispath(cdna_fa_gz_kallisto_index)

        print_and_run_cmd(`kallisto index --index $cdna_fa_gz_kallisto_index $cdna_fa_gz`)

    end

    virus_cdna_fa_gz_kallisto_index = "$virus_cdna_fa_gz.kallisto_index"

    if remake || !ispath(virus_cdna_fa_gz_kallisto_index)

        print_and_run_cmd(`kallisto index --index $virus_cdna_fa_gz_kallisto_index $virus_cdna_fa_gz`)

    end

    Dict(
        "dna.fa.bgz"=>dna_fa_bgz,
        "dna.fa.gz.mmi"=>dna_fa_gz_mmi,
        "chromosome.bed.gz"=>chromosome_bed_gz,
        "chrn_n.tsv"=>chrn_n_tsv,
        "cdna.fa.gz.kallisto_index"=>cdna_fa_gz_kallisto_index,
        "virus_cdna.fa.gz.kallisto_index"=>virus_cdna_fa_gz_kallisto_index,
    )

end
