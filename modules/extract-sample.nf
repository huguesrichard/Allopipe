process EXTRACT_SAMPLE {
	label 'allopipe'
	tag "$sample_id"
	publishDir "${output_dir}/runs/${run_name}/vcf_indiv", mode: 'copy', overwrite: false

	input:
	tuple val(sample_id), path(multi_vcf)
	val  run_name
	val  output_dir

	output:
	tuple val(sample_id),
		path("${sample_id}.vcf.gz"),
		path("${sample_id}.vcf.gz.tbi"),  emit: sample_vcf

	script:
	"""
	bcftools view \
		--samples ${sample_id} \
		${multi_vcf} \
	| bcftools view \
		-e 'GT="mis"' \
		-O z \
		-o ${sample_id}.vcf.gz

	bcftools index --tbi ${sample_id}.vcf.gz
	"""
}
