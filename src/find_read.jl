using Dates


function find_reads(sample_dir::String)
    
    start_time = now()
    
    number_of_fastq_files = 0
    number_of_fastq_gz_files = 0
    fastq_files_to_check = []
    
    println("Walking sample directory...\n")

    for (root, dirs, files) in walkdir("$sample_dir")
        
        println("$root\n")
        
        for file in files
            if occursin(".fastq", file)
                number_of_fastq_files += 1
            end
            if occursin("fastq.gz", file)
                number_of_fastq_gz_files += 1
                push!(fastq_files_to_check, joinpath(root, file))
            end
            if occursin("fq.gz", file)
                number_of_fastq_gz_files += 1
                push!(fastq_files_to_check, joinpath(root, file))
            end
        end
    end

    println("\nNumber of fastq files found in directories walked: $number_of_fastq_files\n")

    println("Number of fastq.gz or fq.gz files found in directories walked: $number_of_fastq_gz_files\n")

    println(string("Number of fastq.gz or fq.gz files to be checked: ", length(fastq_files_to_check)))
    
    end_time = now()
    
    println("\nDone at: $end_time\n")
    
    println("Took $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).\n")
    
    return fastq_files_to_check
    
end


