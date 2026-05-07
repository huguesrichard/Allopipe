process EXTRACT_SAMPLE {
	label 'allopipe'
	tag "$sample_id"

	input:
	tuple val(sample_id), path(multi_vcf)

	output:
	tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: sample_vcf

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
