function check_environment()

    println("Checking environment...")

    for pr::String in (
        "skewer",
        "fastqc",
        "bgzip",
        "tabix",
        "minimap2",
        "samtools",
        "bcftools",
        "kallisto",
    )
        run_command(`which $pr`)

    end

    for pr::String in (
        "configManta.py",
        "configureStrelkaGermlineWorkflow.py",
        "configureStrelkaSomaticWorkflow.py",
    )
        run_command(`bash -c "source activate py2 && which $pr"`)

    end

end

export check_environment