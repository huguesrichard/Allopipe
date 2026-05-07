process FINALIZE_COHORT {
	label 'allopipe'
	publishDir "${output_dir}", mode: 'copy', overwrite: true

	input:
	path run_dirs, stageAs: "pair_runs/??/*"
	val  run_name
	val  output_dir

	output:
	path "runs/${run_name}", emit: results_dir

	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/finalize_cohort.py \
		--run-name ${run_name} \
		--output-dir ./ \
		--run-dir pair_runs/*/*
	"""
}
