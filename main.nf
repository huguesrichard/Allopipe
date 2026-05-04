include { VEP } from './modules/vep-annotation'
include { ALLOCOUNT } from './modules/allo-count'
include { ALLOAFFINITY } from './modules/allo-affinity'

log.info """
    donor         : ${params.donor}
    recipient     : ${params.recipient}
    run_name      : ${params.run_name}
    orientation   : ${params.orientation}
    imputation    : ${params.imputation}
    optional_args : ${params.allocount_opts}
    output_dir    : ${params.output_dir}
""".stripIndent()

workflow AlloPipe {

	take:
	donor_file
	recipient_file

	main:
	// Make VEP annotation optional (--skip_vep_annotation)
	if (!params.skip_vep_annotation) {
		VEP(
			donor_file,
			recipient_file
		)

		// VEP-annotated files as input for ALLOCOUNT
		donor_input     = VEP.out.donor_vep
		recipient_input = VEP.out.recipient_vep

	} else {
		donor_input     = donor_file
		recipient_input = recipient_file
	}


	// Required arguments
	def required = ["donor", "recipient", "run_name", "orientation", "imputation"]
	required.each { p ->
    if (!params[p]) error "ERROR: --${p} is required"
	}

	ALLOCOUNT(
		donor_input,
		recipient_input,
		params.run_name,
		params.orientation,
		params.imputation,
		params.allocount_opts,
		params.output_dir,
	)

	ALLOAFFINITY(
		ALLOCOUNT.out.results_dir, 
		params.run_name,
		params.ensembl_path,
		params.hla_typing,
		params.alloaffinity_opts,
		params.output_dir,
	)	
}

workflow {
    AlloPipe(
        Channel.fromPath(params.donor),
        Channel.fromPath(params.recipient)
    )
}
