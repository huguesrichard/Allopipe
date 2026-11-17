process CHECK_ARGUMENTS_AND_TOOLS {
	label 'allopipe'
	cache false
	debug true

	input:
	val request

	output:
	path 'validated-options.json', emit: options

	script:
	def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
	"""
	PYTHONPATH=${quote("${projectDir}/src")} python -m tools.validate_arguments ${quote(request)}
	PYTHONPATH=${quote("${projectDir}/src")} python -m tools.check_executables validated-options.json
	"""
}
