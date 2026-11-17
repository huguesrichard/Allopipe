process PREPARE_ENSEMBL {
	label 'allopipe'
	container null
	cache false
	tag { "${assembly}, Ensembl ${release}" }
	debug true

	input:
	val assembly
	val release
	val reference_dir
	val vep_version
	val input_vcfs
	val arguments_ready

	output:
	path 'ensembl.ready', emit: ready

	script:
	def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
	def version_arg = vep_version ? "--vep-version ${quote(vep_version)}" : ''
	def vcf_args = input_vcfs.collect { path -> "--vcf ${quote(path)}" }.join(' ')
	"""
	python ${quote("${projectDir}/src/tools/prepare_ensembl.py")} \\
		--assembly ${quote(assembly)} --release ${quote(release)} \\
		--directory ${quote(reference_dir)} ${version_arg} ${vcf_args} 2>&1
	touch ensembl.ready
	"""
}
