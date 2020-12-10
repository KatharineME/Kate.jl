using Dates


function concatenate_reads(
        read_file_paths, 
        sample_name::String, 
        input_dir::String,
    )

    start_time = now()

    forward_read_files = []
    
    reverse_read_files = []

    number_of_forward_reads = 0
    
    number_of_reverse_reads = 0

    for file in read_file_paths
        
        if occursin("R1", file)
            
            push!(forward_read_files, file)
            
            number_of_forward_reads += 1
            
        end
        
    end

    for file in read_file_paths
        
        if occursin("R2", file)
            
            push!(reverse_read_files, file)
            
            number_of_reverse_reads += 1
        end

    end

    println("Number of forward (R1) read files = $number_of_forward_reads\n")
    
    println("Number of reverse (R2) read files = $number_of_reverse_reads\n")

    sample_cat_dir = joinpath(input_dir, string(sample_name, "_cat"))
    
    if ispath(sample_cat_dir)
        
        println("Skipping concatenation because concatenate directory already exists:\n $sample_cat_dir\n")
        
    else

        run(pipeline(`mkdir $sample_cat_dir`))
    
        if number_of_forward_reads > 1

            println("\nCombining R1 reads\n")

            run(pipeline(`cat $forward_read_files`, stdout=joinpath(sample_cat_dir, string(sample_name, "_R1.fastq.gz"))))

        if number_of_reverse_reads > 1

            println("\nCombining R2 reads\n")

            run(pipeline(`cat $reverse_read_files`, stdout=joinpath(sample_cat_dir, string(sample_name, "_R2.fastq.gz"))))       

        else

            println("Number of forward and reverse read files are not more than 1, so there are no fastq files to concatenate.")

        end

    end
        
    end
        
    end_time = now()
            
    println("\nDone at: $end_time\n")
    
    println("Took $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).\n")     
            
end