process ALLO_AFFINITY {
	label 'allopipe'
	tag "$pair_id"
	stageInMode 'copy'
	publishDir "${output_dir}", mode: 'copy', overwrite: true, enabled: params.mode == 'pair'

	input:
	tuple val(pair_id), val(run_name), path(allo_count_results)
	val  ensembl_path
	val  hla_typing
	val  allo_affinity_opts
	val  output_dir

	output:
    tuple val(pair_id), val(run_name), path("runs/${run_name}"), emit: results_dir
	
	script:
	"""
	mkdir -p runs
	cp -R ${allo_count_results} runs/${run_name}

	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/aams_pipeline.py \
		-n ${run_name} \
		-p ${pair_id} \
		-d ${projectDir}/${ensembl_path} \
		-a ${hla_typing} \
		${allo_affinity_opts} \
		-o ./
	"""
}
