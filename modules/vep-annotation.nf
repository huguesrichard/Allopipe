process VEP_ANNOTATION {
	label 'allopipe'
	tag { params.mode == 'cohort' ? sample_id : "${sample_id.toUpperCase()}: ${sample_file.simpleName}" }
	stageInMode 'copy'
	container "ensemblorg/ensembl-vep:${params.vep_version}"

	input:
	tuple val(sample_id), path(sample_file)
	
    output:
    tuple val(sample_id), path("${sample_file.simpleName}_VEP.vcf.gz"), emit: annotated_vcf
	
	script:
	"""
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
		--dir_cache ${params.vep_cache}
	"""
}
