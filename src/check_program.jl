include("print_and_run_cmd.jl")


function check_program()

    println("Checking program ...")

    for program in (
        "skewer",
        "fastqc",
        "bgzip",
        "tabix",
        "minimap2",
        "samtools",
        "configManta.py",
        "configureStrelkaGermlineWorkflow.py",
        "configureStrelkaSomaticWorkflow.py",
        "bcftools",
        "snpEff",
        "kallisto",
    )

        print_and_run_cmd(`which $program`)

    end

end
