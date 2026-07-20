include { VEP_ANNOTATION } from './modules/vep-annotation'
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
		requireParams(['donor', 'recipient', 'hla_typing'])
		return [[pair_id: '', donor: 'donor', recipient: 'recipient', hla_typing: params.hla_typing]]
	}

	requireParams(['multi_vcf', 'pairs'])
	def lines = file(params.pairs).readLines().findAll { it.trim() }
	if (lines.size() < 2) {
		error "ERROR: --pairs must contain a header and at least one donor/recipient row"
	}

	def header = lines[0].split(',', -1).collect { it.trim() }
	def donorIdx = header.indexOf('donor')
	def recipientIdx = header.indexOf('recipient')
	def hlaIdx = header.indexOf('hla')
	if (donorIdx < 0 || recipientIdx < 0) {
		error "ERROR: --pairs must contain donor and recipient columns"
	}
	if (hlaIdx < 0) {
		error "ERROR: --pairs must contain an hla column"
	}
	if (hlaIdx != header.size() - 1) {
		error "ERROR: --pairs must place hla as the last column"
	}

	def body = lines.drop(1)
	def width = body.size().toString().size()
	return body.withIndex().collect { line, index ->
		def cols = line.split(',', 3).collect { it.trim().replaceAll(/^"|"$/, '') }
		if (cols.size() <= Math.max(donorIdx, recipientIdx)) {
			error "ERROR: --pairs row ${index + 2} is missing donor/recipient values"
		}
		def hlaTyping = cols.size() > hlaIdx ? cols[hlaIdx] : null
		if (!hlaTyping?.trim()) {
			error "ERROR: --pairs row ${index + 2} has an empty hla value"
		}
		[
			pair_id: "PAIR_${String.format('%0' + width + 'd', index + 1)}",
			donor: cols[donorIdx],
			recipient: cols[recipientIdx],
			hla_typing: hlaTyping,
		]
	}
}


workflow AlloPipe {
	main:
	if (!(params.mode in ['pair', 'cohort'])) {
		error "ERROR: --mode must be either 'pair' or 'cohort'"
	}
	if (params.mode == 'pair') {
		requireParams(['run_name', 'orientation', 'imputation', 'ensembl_path', 'hla_typing'])
	} else {
		requireParams(['run_name', 'orientation', 'imputation', 'ensembl_path'])
	}
	def frameshiftPlugin = params.frameshift && params.frameshift_plugin_path ?
		file(params.frameshift_plugin_path, checkIfExists: true) :
		file("${projectDir}/modules/vep-annotation.nf")

	def pairRows = buildPairRows()
	pairs_ch = Channel.from(pairRows).map { row -> tuple(row.pair_id, row.donor, row.recipient, row.hla_typing) }

	if (params.mode == 'pair') {
		raw_samples_ch = Channel.of(
			tuple('donor', file(params.donor)),
			tuple('recipient', file(params.recipient)),
		)
		if (!params.skip_vep_annotation) {
			VEP_ANNOTATION(raw_samples_ch.map { sample_id, sample_file -> tuple(sample_id, sample_file, frameshiftPlugin) })
			samples_ch = VEP_ANNOTATION.out.annotated_vcf
		} else {
			samples_ch = raw_samples_ch
		}
	} else {
		def sampleRows = pairRows
			.collectMany { row -> [row.donor, row.recipient] }
			.unique()
			.sort()
			.collect { sample -> tuple(sample, file(params.multi_vcf)) }
		cohort_samples_ch = Channel.from(sampleRows)

		if (!params.skip_vep_annotation) {
			VEP_ANNOTATION(Channel.of(tuple(file(params.multi_vcf).simpleName, file(params.multi_vcf), frameshiftPlugin)))
			def annotated_multi_vcf_ch = VEP_ANNOTATION.out.annotated_vcf.map { cohort_id, annotated_vcf -> annotated_vcf }
			EXTRACT_SAMPLE(
				cohort_samples_ch.combine(annotated_multi_vcf_ch).map { sample_id, multi_vcf, annotated_multi_vcf ->
					tuple(sample_id, annotated_multi_vcf)
				},
				params.run_name,
				params.output_dir,
			)
		} else {
			EXTRACT_SAMPLE(cohort_samples_ch, params.run_name, params.output_dir)
		}
		raw_samples_ch = EXTRACT_SAMPLE.out.sample_vcf.map { sample_id, sample_vcf -> tuple(sample_id, sample_vcf) }
		samples_ch = raw_samples_ch
	}

	donor_join_ch = pairs_ch
		.map { pair_id, donor, recipient, hla_typing -> tuple(donor, pair_id, donor, recipient, hla_typing) }
		.combine(samples_ch, by: 0)
	recipient_join_ch = donor_join_ch
		.map { donor_key, pair_id, donor, recipient, hla_typing, donor_vcf -> tuple(recipient, pair_id, donor, recipient, hla_typing, donor_vcf) }
		.combine(samples_ch, by: 0)
	pair_data_ch = recipient_join_ch
		.map { recipient_key, pair_id, donor, recipient, hla_typing, donor_vcf, recipient_vcf -> tuple(pair_id, donor_vcf, recipient_vcf, hla_typing) }
	pair_inputs_ch = pair_data_ch.map { pair_id, donor_vcf, recipient_vcf, hla_typing -> tuple(pair_id, donor_vcf, recipient_vcf) }

	ALLO_COUNT(
		pair_inputs_ch,
		params.run_name,
		params.orientation,
		params.imputation,
		params.frameshift,
		params.allo_count_opts,
		params.output_dir,
	)

	affinity_inputs_ch = pair_data_ch.join(ALLO_COUNT.out.results_dir, by: 0)

	ALLO_AFFINITY(
		affinity_inputs_ch,
		params.ensembl_path,
		params.allo_affinity_opts,
		params.output_dir,
	)

	if (params.mode == 'cohort') {
		NORMALIZE_COHORT(
			ALLO_AFFINITY.out.results_dir.map { pair_id, run_name, run_dir -> run_dir }.collect(),
			EXTRACT_SAMPLE.out.sample_vcf.flatMap { sample_id, sample_vcf -> [sample_vcf] }.collect(),
			params.run_name,
			params.output_dir,
		)
	}
}


workflow {
	AlloPipe()
}
