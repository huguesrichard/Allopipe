process VEP_ANNOTATION {
	label 'allopipe'
	tag { params.mode == 'cohort' ? sample_id : "${sample_id.toUpperCase()}: ${sample_file.simpleName}" }
	stageInMode 'copy'
	container "ensemblorg/ensembl-vep:${params.vep_version}"

	input:
	tuple val(sample_id), path(sample_file), path(frameshift_plugin_path)
	
    output:
    tuple val(sample_id), path("${sample_file.simpleName}_VEP.vcf.gz"), emit: annotated_vcf
	
	script:
	def frameshift_plugin = params.frameshift ? '--plugin Frameshift' : ''
	def frameshift_plugin_dir = ''
	if (params.frameshift && params.frameshift_plugin_path) {
		def dir_plugins = frameshift_plugin_path.getName().endsWith('.pm') ? '.' : frameshift_plugin_path
		frameshift_plugin_dir = "--dir_plugins ${dir_plugins}"
	}
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
		${frameshift_plugin} \
		${frameshift_plugin_dir} \
		--dir_cache ${params.vep_cache}
	"""
}
