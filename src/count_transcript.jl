using Dates

include("print_and_run_cmd.jl")


function count_transcript(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    cdna_fasta_gz::String,
    output_dir::String,
    n_job::Int,
)

    start_time = now()

    println("($start_time) Counting transcript ...")

    cdna_fasta_gz_kallisto_index::String = "$cdna_fasta_gz.kallisto_index"

    if !ispath(cdna_fasta_gz_kallisto_index)

        print_and_run_cmd(`kallisto index --index $cdna_fasta_gz_kallisto_index $cdna_fasta_gz`)

    end

    mkpath(output_dir)

    print_and_run_cmd(`kallisto quant --threads $n_job --index $cdna_fasta_gz_kallisto_index --output-dir $output_dir $_1_fastq_gz $_2_fastq_gz`)

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end
