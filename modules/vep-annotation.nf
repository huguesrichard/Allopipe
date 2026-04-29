process VEP {
	label 'allopipe'
	stageInMode 'copy'

	input:
	path donor_file
	path recipient_file
	

    output:
    path "${donor_file.simpleName}_VEP.vcf.gz", emit: donor_vep
    path "${recipient_file.simpleName}_VEP.vcf.gz", emit: recipient_vep

	
	script:
	"""	
	module load singularity
	
	# Donor
	singularity exec ${params.sif_dir}/vep.sif \
    vep --dir ${HOME}/vep_data --cache --assembly GRCh38 --offline --af_gnomade \
        --fork 4 --coding_only --pick_allele --use_given_ref --vcf --compress_output gzip \
   	 	-i ${donor_file} \
		-o ${donor_file.simpleName}_VEP.vcf.gz

	# Recipient
	singularity exec ${params.sif_dir}/vep.sif \
    vep --dir $HOME/vep_data --cache --assembly GRCh38 --offline --af_gnomade \
        --fork 4 --coding_only --pick_allele --use_given_ref --vcf --compress_output gzip \
		-i ${recipient_file} \
		-o ${recipient_file.simpleName}_VEP.vcf.gz
	"""
}