include("print_command_and_run.jl")

function make_grch38_dna_fasta_for_alignment(
    directory_path::String;
    overwrite=false,
)::String

    data_grch_directory_path = "$DATA_DIRECTORY_PATH/bwa.kit/resource-GRCh38"

    final_fa_file_path = "$directory_path/GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa"

    final_fa_gz_file_path = "$final_fa_file_path.gz"

    if !overwrite && isfile(final_fa_gz_file_path):

        error(final_fa_gz_file_path)

    fa_gz_file_path = "$directory_path/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"

    download(
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/$(split(fa_gz_file_path, "/")[end])",
        directory_path,
    )

    print_command_and_run(pipeline(
        `gzip --decompress --to-stdout $fa_gz_file_path`,
        stdout=final_fa_file_path,
    ))

    print_command_and_run(pipeline(
        `cat $data_grch_directory_path/hs38DH-extra.fa`,
        stdout=final_fa_file_path,
        append=true,
    ))

    print_command_and_run(pipeline(
        `gzip $final_fa_file_path --to-stdout`,
        stdout=final_fa_gz_file_path,
    ))

    for file_path = (
        final_fa_file_path,
        final_fa_gz_file_path,
    )

        print_command_and_run(`cp $data_grch_directory_path/hs38DH.fa.alt $file_path.alt`)

    return final_fa_gz_file_path

end
