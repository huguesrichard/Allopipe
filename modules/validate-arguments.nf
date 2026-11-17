process VALIDATE_ARGUMENTS {
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
	"""
}
