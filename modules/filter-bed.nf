process FILTER_BED {
	label 'allopipe'
	tag "$sample_id"
	publishDir "${output_dir}/runs/${run_name}/vcf_vep", mode: 'copy', overwrite: true,
		saveAs: { filename -> params.mode == 'pair' ? "${sample_id}_${filename}" : filename }

	input:
	tuple val(sample_id), path(sample_file)
	path bed_file
	val run_name
	val output_dir

	output:
	tuple val(sample_id), path("${sample_file.simpleName}_bed.vcf.gz"), emit: filtered_vcf
	path "${sample_file.simpleName}_bed.vcf.gz.tbi", emit: index

	script:
	"""
	bcftools view --targets-file "${bed_file}" --targets-overlap 1 \\
		-O z -o "${sample_file.simpleName}_bed.vcf.gz" "${sample_file}"
	bcftools index --tbi "${sample_file.simpleName}_bed.vcf.gz"
	"""
}
