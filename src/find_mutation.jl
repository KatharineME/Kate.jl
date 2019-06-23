function find_mutation(
        germ_bam::String,
        soma_bam::String,
        dna_fa_bgz::String,
        chromosome_bed_gz::String,
        sequencing_scope::String,
        n_job::Int,
        gb_memory::Int,
        manta_dir::String,
        strelka_dir::String,
        snpeff_dir::String,
    )

    println("Finding mutation ...")

    config_parameters = "--referencefasta $dna_fa_bgz --callregions $chromosome_bed_gz --$sequencing_scope"

    if ispath(germ_bam) && ispath(soma_bam)

        config_parameters = "$config_parameters --normalbam $germ_bam --tumorbam $soma_bam"

    elseif ispath(germ_bam)

        config_parameters = "$config_parameters --bam $germ_bam"

    else

        error("Arguments do not contain either germ and soma .bams or germ .bam.")

    end

  run_parameters = "--mode local --jobs $n_job --memGb $gb_memory --quiet"

  Kraft.print_and_run_cmd(`source activate py2.7`)

  Kraft.print_and_run_cmd(`configManta.py $config_parameters --outputContig --runDir $manta_dir`)

  Kraft.print_and_run_cmd(`$manta_dir/runWorkflow.py $run_parameters`)

  if [[ -e $GERM_BAM && -e $SOMA_BAM ]]; then

    configureStrelkaSomaticWorkflow.py $CONFIG_PARAMETERS --runDir $STRELKA_DIR --indelCandidates $MANTA_DIR/results/variants/candidateSmallIndels.vcf.gz

  else

    configureStrelkaGermlineWorkflow.py $CONFIG_PARAMETERS --runDir $STRELKA_DIR

  fi

  $STRELKA_DIR/runWorkflow.py $RUN_PARAMETERS

  conda deactivate

  if [[ -e $GERM_BAM && -e $SOMA_BAM ]]; then

    local INDEL_VCF_GZ=$STRELKA_DIR/results/variants/somatic.indels.vcf.gz

    local SNV_VCF_GZ=$STRELKA_DIR/results/variants/somatic.snvs.vcf.gz

    local STRELKA_SOMA_SAMPLE_NAME_TXT=/tmp/strelka_soma_sample_name.txt

    touch $STRELKA_SOMA_SAMPLE_NAME_TXT

    echo "Germ" >> $STRELKA_SOMA_SAMPLE_NAME_TXT

    echo "Soma" >> $STRELKA_SOMA_SAMPLE_NAME_TXT

    bcftools reheader --threads $N_JOB --samples $STRELKA_SOMA_SAMPLE_NAME_TXT $INDEL_VCF_GZ > /tmp/indel.vcf.gz

    mv --force /tmp/indel.vcf.gz $INDEL_VCF_GZ

    tabix --force $INDEL_VCF_GZ

    bcftools reheader --threads $N_JOB --samples $STRELKA_SOMA_SAMPLE_NAME_TXT $SNV_VCF_GZ > /tmp/snv.vcf.gz

    mv --force /tmp/snv.vcf.gz $SNV_VCF_GZ

    tabix --force $SNV_VCF_GZ

    rm --force $STRELKA_SOMA_SAMPLE_NAME_TXT

    local CONCAT_VCFS="$MANTA_DIR/results/variants/somaticSV.vcf.gz $INDEL_VCF_GZ $SNV_VCF_GZ"

  else

    local CONCAT_VCFS="$MANTA_DIR/results/variants/diploidSV.vcf.gz $STRELKA_DIR/results/variants/variants.vcf.gz"

  fi

  bcftools concat --allow-overlaps --threads $N_JOB $CONCAT_VCFS | bcftools annotate --rename-chrs $CHRN_N_TSV --threads $N_JOB | bgzip --stdout --threads $N_JOB > /tmp/concat.vcf.gz

  tabix /tmp/concat.vcf.gz

  snpEff -Xms${GB_MEMORY}g -Xmx${GB_MEMORY}g GRCh38.86 -verbose -noLog -csvStats $SNPEFF_DIR/stats.csv -htmlStats $SNPEFF_DIR/stats.html /tmp/concat.vcf.gz | bgzip --stdout --threads $N_JOB > $SNPEFF_DIR/variant.vcf.gz

  tabix $SNPEFF_DIR/variant.vcf.gz

  rm --force /tmp/concat.vcf.gz /tmp/concat.vcf.gz.tbi

end
