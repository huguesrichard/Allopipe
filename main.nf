include { VALIDATE_ARGUMENTS } from './modules/validate-arguments'
include { PREPARE_ENSEMBL } from './modules/prepare-ensembl'
include { FILTER_BED } from './modules/filter-bed'
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

def requireVepAnnotation(path) {
	def stream = java.nio.file.Files.newInputStream(path)
	try {
		def input = path.toString().endsWith('.gz') ? new java.util.zip.GZIPInputStream(stream) : stream
		def reader = new BufferedReader(new InputStreamReader(input, 'UTF-8'))
		try {
			String line
			while ((line = reader.readLine()) != null && !line.startsWith('#CHROM')) {
				if (line.startsWith('##VEP-command-line=')) return
			}
		} finally {
			reader.close()
		}
	} finally {
		stream.close()
	}
	error "ERROR: --skip_vep_annotation requires a VEP-annotated VCF; missing ##VEP-command-line header in ${path}"
}

def checkPreRun(inputVcfs) {
	inputVcfs.each { requireVepAnnotation(file(it)) }
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
			pair_id: "P${String.format('%0' + width + 'd', index + 1)}",
			donor: cols[donorIdx],
			recipient: cols[recipientIdx],
			hla_typing: hlaTyping,
		]
	}
}

def checkedOutputRunPath(runName, outputDir) {
	def runsPath = file("${outputDir}/runs").toAbsolutePath().normalize()
	def outputRunPath = runsPath.resolve(runName.toString()).normalize()

	if (outputRunPath.parent != runsPath) {
		error "ERROR: --run_name must be a single directory name, not a path:\n${runName}\n\n" +
			"AlloPipe automatically creates the run directory at:\n" +
			"${runsPath}/<run_name>\n\n" +
			"Use --output_dir to choose a different parent output directory."
	}
	return outputRunPath
}

def prepareOutputRun(runName, outputDir, forceOverwrite, resumeRun) {
	def outputRunPath = checkedOutputRunPath(runName, outputDir)
	if (!outputRunPath.exists()) {
		return
	}
	if (resumeRun && !forceOverwrite) {
		log.info "RESUME: keeping existing output run directory:\n${outputRunPath}"
		return
	}
	if (!forceOverwrite) {
		error "ERROR: output run directory already exists:\n${outputRunPath}\n\n" +
			"Solutions:\n" +
			"  - Replace it: pass --force_overwrite\n" +
			"  - Continue it: use Nextflow -resume\n" +
			"  - Keep it: choose a unique combination of --output_dir and --run_name"
	}
	if (java.nio.file.Files.isSymbolicLink(outputRunPath) || !outputRunPath.toFile().isDirectory()) {
		error "ERROR: refusing to overwrite output run path because it is not a regular directory:\n${outputRunPath}"
	}

	log.warn "--force_overwrite is removing existing output run directory:\n${outputRunPath}"
	if (!outputRunPath.toFile().deleteDir() || outputRunPath.exists()) {
		error "ERROR: failed to remove existing output run directory:\n${outputRunPath}"
	}
}


workflow AlloPipe {
	main:
	if (!(params.mode in ['pair', 'cohort'])) {
		error "ERROR: --mode must be either 'pair' or 'cohort'"
	}
	if (params.mode == 'pair') {
		requireParams(['run_name', 'orientation', 'imputation', 'assembly', 'hla_typing'])
		} else {
			requireParams(['run_name', 'orientation', 'imputation', 'assembly'])
			if (params.hla_typing) {
				log.warn "--hla_typing is ignored in cohort mode; per-pair HLA typing is read from the mandatory 'hla' column in the CSV file provided with --pairs"
			}
		}
	if (!(params.assembly in ['GRCh37', 'GRCh38'])) {
		error "ERROR: --assembly must be GRCh37 or GRCh38"
	}
	def ensemblRelease = params.ensembl_release
	if (ensemblRelease == null) {
		def vepTag = params.vep_version?.toString() ?: ''
		def match = vepTag =~ /^release_(\d+)(?:\.\d+)?$/
		if (!match.matches()) {
			error "ERROR: cannot derive Ensembl release from --vep_version '${vepTag}'; specify --ensembl_release"
		}
		ensemblRelease = match.group(1)
	}
	if (!(ensemblRelease.toString() ==~ /[0-9]+/) || ensemblRelease.toBigInteger() < 1) {
		error "ERROR: --ensembl_release must be a positive integer"
	}
	def publishedOutputDir = file(params.output_dir).toAbsolutePath().normalize().toString()
	checkedOutputRunPath(params.run_name, publishedOutputDir)
	def nextflowCommandBase64 = workflow.commandLine.bytes.encodeBase64().toString()
	def bedFile = params.bed_filter ? file(params.bed_file ?: "${projectDir}/tutorial/Twist_Exome_Core_Covered_Targets_hg38.bed", checkIfExists: true) : null
	def frameshiftPlugin = params.frameshift && params.frameshift_plugin_path ?
		file(params.frameshift_plugin_path, checkIfExists: true) :
		file("${projectDir}/modules/vep-annotation.nf")
	def vepAssembly = params.assembly
	def ensemblPath = "${projectDir}/data/Ensembl/${vepAssembly}".toString()

	def pairRows = buildPairRows()
	def inputVcfs = params.skip_vep_annotation ?
		(params.mode == 'pair' ? [params.donor, params.recipient] : [params.multi_vcf])
			.collect { file(it, checkIfExists: true).toAbsolutePath().toString() } : []
	checkPreRun(inputVcfs)
	def validationRequest = groovy.json.JsonOutput.toJson([
		run_name: params.run_name, orientation: params.orientation, imputation: params.imputation,
		mode: params.mode, frameshift: params.frameshift, ensembl_path: ensemblPath,
		count_opts: params.allo_count_opts ?: '', affinity_opts: params.allo_affinity_opts ?: '',
		pairs: pairRows.collect { [pair_id: it.pair_id, hla_typing: it.hla_typing] },
	]).bytes.encodeBase64().toString()
	VALIDATE_ARGUMENTS(validationRequest)
	def validatedArguments = VALIDATE_ARGUMENTS.out.options.map { optionsFile ->
		prepareOutputRun(params.run_name, publishedOutputDir, params.force_overwrite, workflow.resume)
		def options = new groovy.json.JsonSlurper().parseText(optionsFile.text)
		tuple(optionsFile, options.count_opts, options.affinity_opts)
	}
	def argumentsReady = validatedArguments.map { marker, countOptions, affinityOptions -> marker }
	def validatedCountOptions = validatedArguments.map { marker, countOptions, affinityOptions -> countOptions }
	def validatedAffinityOptions = validatedArguments.map { marker, countOptions, affinityOptions -> affinityOptions }
	PREPARE_ENSEMBL(
		params.assembly,
		ensemblRelease,
		ensemblPath,
		params.skip_vep_annotation ? '' : params.vep_version,
		inputVcfs,
		argumentsReady,
	)
	def referencesReady = PREPARE_ENSEMBL.out.ready
	pairs_ch = Channel.from(pairRows).map { row -> tuple(row.pair_id, row.donor, row.recipient, row.hla_typing) }

	if (params.mode == 'pair') {
		raw_samples_ch = Channel.of(
			tuple('donor', file(params.donor)),
			tuple('recipient', file(params.recipient)),
		).combine(referencesReady).map { sample_id, sample_file, ready -> tuple(sample_id, sample_file) }
		if (params.bed_filter) {
			FILTER_BED(raw_samples_ch, bedFile, params.run_name, publishedOutputDir)
			raw_samples_ch = FILTER_BED.out.filtered_vcf
		}
		if (!params.skip_vep_annotation) {
			VEP_ANNOTATION(
				raw_samples_ch.map { sample_id, sample_file -> tuple(sample_id, sample_file, frameshiftPlugin, vepAssembly) },
				params.run_name,
				publishedOutputDir,
			)
			samples_ch = VEP_ANNOTATION.out.annotated_vcf
		} else {
			samples_ch = raw_samples_ch
		}
	} else {
		cohort_input_ch = Channel.of(tuple(file(params.multi_vcf).simpleName, file(params.multi_vcf)))
			.combine(referencesReady).map { cohort_id, cohort_vcf, ready -> tuple(cohort_id, cohort_vcf) }
		if (params.bed_filter) {
			FILTER_BED(cohort_input_ch, bedFile, params.run_name, publishedOutputDir)
			cohort_input_ch = FILTER_BED.out.filtered_vcf
			cohort_bed_vcfs = FILTER_BED.out.filtered_vcf.map { cohort_id, vcf -> vcf }.mix(FILTER_BED.out.index).collect()
		} else {
			cohort_bed_vcfs = Channel.value([])
		}
		def sampleRows = pairRows
			.collectMany { row -> [row.donor, row.recipient] }
			.unique()
			.sort()
			.collect { sample -> sample }
		cohort_samples_ch = Channel.from(sampleRows)

		if (!params.skip_vep_annotation) {
			VEP_ANNOTATION(
				cohort_input_ch.map { cohort_id, cohort_vcf -> tuple(cohort_id, cohort_vcf, frameshiftPlugin, vepAssembly) },
				params.run_name,
				publishedOutputDir,
			)
			def annotated_multi_vcf_ch = VEP_ANNOTATION.out.annotated_vcf.map { cohort_id, annotated_vcf -> annotated_vcf }
			cohort_vep_vcfs = annotated_multi_vcf_ch.mix(VEP_ANNOTATION.out.report).collect()
			EXTRACT_SAMPLE(
				cohort_samples_ch.combine(annotated_multi_vcf_ch).map { sample_id, annotated_multi_vcf ->
					tuple(sample_id, annotated_multi_vcf)
				},
				params.run_name,
				publishedOutputDir,
			)
		} else {
			cohort_vep_vcfs = Channel.value([])
			EXTRACT_SAMPLE(cohort_samples_ch.combine(cohort_input_ch.map { cohort_id, cohort_vcf -> cohort_vcf }), params.run_name, publishedOutputDir)
		}
		raw_samples_ch = EXTRACT_SAMPLE.out.sample_vcf.map { sample_id, sample_vcf, sample_tbi -> tuple(sample_id, sample_vcf) }
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
		validatedCountOptions,
		publishedOutputDir,
		nextflowCommandBase64,
	)

	affinity_inputs_ch = pair_data_ch.join(ALLO_COUNT.out.results_dir, by: 0)

	ALLO_AFFINITY(
		affinity_inputs_ch,
		ensemblPath,
		validatedAffinityOptions,
		publishedOutputDir,
	)

	if (params.mode == 'cohort') {
		NORMALIZE_COHORT(
			ALLO_AFFINITY.out.results_dir.map { pair_id, run_name, run_dir -> run_dir }.collect(),
			EXTRACT_SAMPLE.out.sample_vcf.flatMap { sample_id, sample_vcf, sample_tbi -> [sample_vcf, sample_tbi] }.collect(),
			cohort_vep_vcfs,
			cohort_bed_vcfs,
			params.run_name,
			publishedOutputDir,
		)
	}
}


workflow {
	AlloPipe()
}
