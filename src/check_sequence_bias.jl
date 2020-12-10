using Dates

function check_sequence_bias(sample_name::String, output_dir::String)
    
    start_time = now()
    
    check_sequence_directory = joinpath(output_dir, string("check_sequence_", sample_name))

    ProcessSequence.print_and_run_cmd(`multiqc --outdir $check_sequence_directory $check_sequence_directory`)
    
    end_time = now()
    
    println("\nDone at: $end_time\n")
    
    println("Took $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).\n")
    
end