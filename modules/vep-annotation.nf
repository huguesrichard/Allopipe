process VEP_ANNOTATION {
	label 'allopipe'
	tag "$sample_id"
	stageInMode 'copy'

	container 'ensemblorg/ensembl-vep:latest'							// specify version (+ cache)

	input:
	tuple val(sample_id), path(sample_file)
	

    output:
    tuple val(sample_id), path("${sample_file.simpleName}_VEP.vcf.gz"), emit: annotated_vcf

	
	script:
	"""
	# Docker (local)
	vep \
		-i ${sample_file} \
		-o ${sample_file.simpleName}_VEP.vcf.gz \
		--vcf \
		--fork 4 \
		--cache \
		--offline \
		--compress_output gzip \
		--force_overwrite \
		--assembly GRCh38 \
		--af_gnomade \
		--coding_only \
		--pick_allele \
		--use_given_ref	\
		--dir_cache /cache	
	"""
	// ###################################################
	// ## Singularity (Slurm)
	// #module load singularity
	
	// #singularity exec ${params.sif_dir}/vep.sif \
	// #vep --dir ${HOME}/vep_data \
	// ###################################################
}
