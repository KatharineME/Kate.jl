function check_programs()

    println("Checking programs ...")

    for program in (
            "skewer",
            "fastqc",
            "minimap2",
            "bgzip",
            "samtools",
            "tabix",
            "bcftools",
            "snpEff",
            "kallisto",
        )

        Kraft.print_and_run_cmd(`which $program`)

    end

end
