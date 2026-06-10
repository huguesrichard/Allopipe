include { VEP } from './modules/vep-annotation'
include { EXTRACT_SAMPLE } from './modules/extract-sample'
include { ALLO_COUNT } from './modules/allo-count'
include { ALLO_AFFINITY } from './modules/allo-affinity'
include { NORMALIZE_COHORT } from './modules/normalize-cohort'


def requireParams(names) {
	names.each { name ->
		if (!params[name]) {
			error "ERROR: --${name} is required"
		}
	}
}


def buildPairRows() {
	if (params.mode == 'pair') {
		requireParams(['donor', 'recipient'])
		return [[pair_id: '', donor: 'donor', recipient: 'recipient']]
	}

	requireParams(['multi_vcf', 'pairs'])
	def lines = file(params.pairs).readLines().findAll { it.trim() }
	if (lines.size() < 2) {
		error "ERROR: --pairs must contain a header and at least one donor/recipient row"
	}

	def header = lines[0].split(',').collect { it.trim() }
	def donorIdx = header.indexOf('donor')
	def recipientIdx = header.indexOf('recipient')
	if (donorIdx < 0 || recipientIdx < 0) {
		error "ERROR: --pairs must contain donor and recipient columns"
	}

	def body = lines.drop(1)
	def width = body.size().toString().size()
	return body.withIndex().collect { line, index ->
		def cols = line.split(',').collect { it.trim() }
		[
			pair_id: "P${String.format('%0' + width + 'd', index + 1)}",
			donor: cols[donorIdx],
			recipient: cols[recipientIdx],
		]
	}
}


workflow AlloPipe {
	main:
	requireParams(['run_name', 'orientation', 'imputation', 'ensembl_path', 'hla_typing'])
	if (!(params.mode in ['pair', 'cohort'])) {
		error "ERROR: --mode must be either 'pair' or 'cohort'"
	}

	def pairRows = buildPairRows()
	pairs_ch = Channel.from(pairRows).map { row -> tuple(row.pair_id, row.donor, row.recipient) }

	if (params.mode == 'pair') {
		raw_samples_ch = Channel.of(
			tuple('donor', file(params.donor)),
			tuple('recipient', file(params.recipient)),
		)
	} else {
		def sampleRows = pairRows
			.collectMany { row -> [row.donor, row.recipient] }
			.unique()
			.sort()
			.collect { sample -> tuple(sample, file(params.multi_vcf)) }
		EXTRACT_SAMPLE(Channel.from(sampleRows), params.run_name, params.output_dir)
		raw_samples_ch = EXTRACT_SAMPLE.out.sample_vcf.map { sample_id, sample_vcf, sample_vcf_index -> tuple(sample_id, sample_vcf) }
	}

	if (!params.skip_vep_annotation) {
		VEP(raw_samples_ch)
		samples_ch = VEP.out.annotated_vcf
	} else {
		samples_ch = raw_samples_ch
	}

	donor_join_ch = pairs_ch
		.map { pair_id, donor, recipient -> tuple(donor, pair_id, donor, recipient) }
		.join(samples_ch)
	recipient_join_ch = donor_join_ch
		.map { donor_key, pair_id, donor, recipient, donor_vcf -> tuple(recipient, pair_id, donor, recipient, donor_vcf) }
		.join(samples_ch)
	pair_inputs_ch = recipient_join_ch
		.map { recipient_key, pair_id, donor, recipient, donor_vcf, recipient_vcf -> tuple(pair_id, donor_vcf, recipient_vcf) }

	ALLO_COUNT(
		pair_inputs_ch,
		params.run_name,
		params.orientation,
		params.imputation,
		params.allo_count_opts,
		params.output_dir,
	)

	ALLO_AFFINITY(
		pair_inputs_ch,
		ALLO_COUNT.out.results_dir,
		params.ensembl_path,
		params.hla_typing,
		params.allo_affinity_opts,
		params.output_dir,
	)

	if (params.mode == 'cohort') {
		NORMALIZE_COHORT(
			ALLO_AFFINITY.out.results_dir.map { pair_id, run_name, run_dir -> run_dir }.collect(),
			EXTRACT_SAMPLE.out.sample_vcf.flatMap { sample_id, sample_vcf, sample_vcf_index -> [sample_vcf, sample_vcf_index] }.collect(),
			params.run_name,
			params.output_dir,
		)
	}
}


workflow {
	AlloPipe()
}
