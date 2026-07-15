process NORMALIZE_COHORT {
	label 'allopipe'
	tag { params.mode == 'cohort' ? "${run_dirs.size()} PAIRS" : null }
	publishDir "${output_dir}", mode: 'copy', overwrite: true

	input:
	path run_dirs, stageAs: "pair_runs/??/*"
	path extracted_vcfs, stageAs: "runs/${params.run_name}/vcf_indiv/*"
	val  run_name
	val  output_dir

	output:
	path "runs/${run_name}", emit: results_dir

	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/tools/normalize_cohort.py \
		--run-name ${run_name} \
		--output-dir ./ \
		--run-dir pair_runs/*/*
	"""
}
