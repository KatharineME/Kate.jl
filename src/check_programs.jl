include("print_and_run_cmd.jl")


function check_programs()

    println("Checking programs ...")

    for program::String in (
        "skewer",
        "fastqc",
        "bgzip",
        "tabix",
        "minimap2",
        "samtools",
        "bcftools",
        "snpEff",
        "kallisto",
    )

        print_and_run_cmd(`which $program`)

    end

    for program::String in (
        "configManta.py",
        "configureStrelkaGermlineWorkflow.py",
        "configureStrelkaSomaticWorkflow.py",
    )

        print_and_run_cmd(`bash -c "source activate py2.7 && which $program"`)

    end

    nothing

end
