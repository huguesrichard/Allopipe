process EXTRACT_SAMPLE {
	label 'allopipe'
	tag "$sample_id"
	publishDir "${output_dir}/runs/${run_name}/vcf_indiv", mode: 'copy', overwrite: true

	input:
	tuple val(sample_id), path(multi_vcf)
	val  run_name
	val  output_dir

	output:
	tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: sample_vcf
	path "${sample_id}.vcf.gz.tbi", emit: sample_vcf_index

	script:
	"""
	bcftools view \
		--samples ${sample_id} \
		--output-type z \
		--output-file ${sample_id}.vcf.gz \
		${multi_vcf}

	bcftools index --tbi ${sample_id}.vcf.gz
	"""
}
