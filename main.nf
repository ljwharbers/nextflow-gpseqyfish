// Include samplesheet plugin
include { fromSamplesheet } from 'plugin/nf-validation'

// Print pipeline info
log.info """\

		G P S E Q   Y F I S H    P I P E L I N E
		========================================
		Indir			: ${params.indir}
		Outdir			: ${params.outdir}
		"""
		.stripIndent()

// Convert .nd2 file to .tiff

process ND2_TO_TIFF {
	label "process_low"
	tag "Converting ${nd2} to .tiff files"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(nd2)
		val dz

	output:
		tuple val(sample), path("${sample}/${params.dapi}*.tiff"), emit: dapi_tiffs
		tuple val(sample), path("${sample}/${params.yfish}*.tiff"), emit: yfish_tiffs
		tuple val(sample), path("${sample}/nd2_to_tiff.log.txt"), emit: log

	script:
		"""
		radiantkit nd2_to_tiff ${nd2} -Z ${dz}
		"""
}

// Find out of focus
process FIND_OOF {
	label "process_high"
	tag "Finding out of focus in DAPI .tiff files from ${sample}"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(dapi_tiffs)

	output:
		tuple val(sample), path("oof.tsv"), emit: tsv
		tuple val(sample), path("oof.args.pkl"), emit: oof_args
		tuple val(sample), path("oof.log.txt"), emit: oof_log
		tuple val(sample), path("*.tiff", includeInputs: true), emit: dapi_infocus

	script:
		"""
		radiantkit tiff_findoof . --rename --threads ${task.cpus}
		"""
}

// Segment TIFF files
process SEGMENT_TIFF {
	label "process_high"
	tag "Segmenting ${sample} tiff files"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(dapi_infocus)

	output:
		tuple val(sample), path("*mask.tiff"), emit: dapi_masks
		tuple val(sample), path("tiff_segment.log.txt"), emit: mask_log

	script:
		"""
		radiantkit tiff_segment .
		--TCZYX \
		--threads ${task.cpus} \
               	--gaussian 2.0 \
              	 --inreg "^${dapi}.*\.tif$" \
               	-y
		"""
}

// Measure objects
process MEASURE_OBJECTS {
	label "process_high"
	tag "Measuring objects in ${sample} tiff files"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(dapi_infocus)
		tuple val(sample), path(dapi_mask)
		val dapi
		val dx
		val dy
		val dz

	output:
		tuple val(sample), path("objects/nuclear_features.tsv")

	script:
		"""
		radiantkit measure_objects . '${dapi}' --aspect ${dz} ${dy} ${dx} --threads ${task.cpus} -y
		"""
}

// Select nuclei
process SELECT_NUCLEI {
	label "process_medium"
	tag "Selecting nuclei in ${tiff}*.tiff"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(dapi_infocus)
		tuple val(sample), path(dapi_mask)
		val k_sigma
		val dapi

	output:
		tuple val(sample), path("*mask_selected.tiff"), emit: mask_selected
		tuple val(sample), path("select_nuclei.log.txt"), emit: select_log
		tuple val(sample), path("select_nuclei.data.tsv"), emit: select_data
		tuple val(sample), path("select_nuclei.args.pkl"), emit: pkl_args
		tuple val(sample), path("select_nuclei.fit.pkl"), emit: pkl_fit


	script:
		"""
		radiantkit select_nuclei . '${dapi}' --k-sigma ${k_sigma} --threads ${task.cpus} -y
		"""
}

// Get radial population
process RADIAL_POPULATION {
	label "process_high"
	tag "Getting radial population of ${sample}"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	input:
		tuple val(sample), path(dapi_infocus)
		tuple val(sample), path(mask_selected)
		tuple val(sample), path(dapi_mask)
		tuple val(sample), path(yfish_tiffs)
		val dapi
		val dx
		val dy
		val dz

	output:
		tuple val(sample), path("radial_population.args.pkl"), emit: radial_pop
		tuple val(sample), path("radial_population.log.txt"), emit: radial_log
		tuple val(sample), path("objects"), emit: objects

	script:
		"""
		radiantkit radial_population . '${dapi}' --aspect ${dz} ${dx} ${dy} --mask-suffix mask_selected --threads ${task.cpus} --slice2d -y
		"""
}

// Create directory structure
process CREATE_DIRS {
	label "process_single"
	tag "Generating directory structure for report"

	input:
		tuple val(sample),
			  path(oof),
			  path(oof_args),
			  path(oof_log),
			  path(select_data),
			  path(select_pkl_args),
			  path(select_pkl_fit),
			  path(select_log),
			  path(radial_pop),
			  path(radial_log),
			  path(radial_objects),
			  path(objects),
			  path(a),
			  path(b),
			  path(c),
			  path(d)

	output:
		path "${sample}/"

	script:
		"""
		mkdir -p ${sample}
		mv -t ${sample}/ ${oof} ${oof_args} ${oof_log} ${select_data} ${select_pkl_args} \
		${select_pkl_fit} ${select_log} ${radial_pop} ${radial_log} ${radial_objects} \
		${objects} ${a} ${b} ${c} ${d}
		"""
}

// Get Radiant report
process RADIANT_REPORT {
	label "process_high"
	tag "Generating Radiant report for all samples"

	container "library://ljwharbers/gpseq/radiantkit:0.0.2"

	publishDir "${params.outdir}", mode: 'copy'

	input:
		path all

	output:
		path "radiant.report.html"

	script:
		"""
		radiantkit report .
		"""
}

workflow {
	// Get .nd2 files available in the input directory
	if (params.tif) {
		Channel
			.fromPath(params.indir + "/**/*.tif*")
			.map { file -> tuple(file.baseName, file) }
			.set { tiff }
		tiff.view()
	} else {
		Channel
			.fromPath(params.indir + "/*.nd2")
			.map { file -> tuple(file.baseName, file) }
			.set { nd2_files }
		tiff = ND2_TO_TIFF(nd2_files, params.dz)
	}

	oof = FIND_OOF(tiff.dapi_tiffs)
	segmented = SEGMENT_TIFF(oof.dapi_infocus)
	objects = MEASURE_OBJECTS(oof.dapi_infocus, segmented.dapi_masks,
							  params.dapi, params.dx, params.dy, params.dz)
	select = SELECT_NUCLEI(oof.dapi_infocus, segmented.dapi_masks, params.k_sigma, params.dapi)
	radial = RADIAL_POPULATION(oof.dapi_infocus,
							   select.mask_selected,
							   segmented.dapi_masks,
							   tiff.yfish_tiffs,
							   params.dapi, params.dx,
							   params.dy, params.dz)

	// Get all output required for report in correct structure
	oof.tsv
		.join(oof.oof_args)
		.join(oof.oof_log)
		.join(select.select_data)
		.join(select.pkl_args)
		.join(select.pkl_fit)
		.join(select.select_log)
		.join(radial.radial_pop)
		.join(radial.radial_log)
		.join(radial.objects)
		.join(objects)
		.join(tiff.dapi_tiffs)
		.join(tiff.yfish_tiffs)
		.join(segmented.dapi_masks)
		.join(select.mask_selected)
		.set { report_input }

	dirs = CREATE_DIRS(report_input)
	RADIANT_REPORT(dirs.collect())
}
