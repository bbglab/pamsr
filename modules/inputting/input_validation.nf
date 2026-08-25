// modules/signature_extraction/input_validation.nf

def getSep(delimiter) {
    def d = delimiter.toLowerCase()
    if (d == 'csv') {
        return ','
    } else if (d == 'tsv') {
        return '\t'
    } else {
        throw new IllegalArgumentException("Unsupported file delimiter format: ${delimiter}")
    }
}

workflow VALIDATE_INPUTS {
    take:
    metadata_path
    metadata_delim
    samplesheet_path
    samplesheet_delim

    main:
    // Parse & validate metadata
    ch_metadata = Channel
        .fromPath(metadata_path, checkIfExists: true)
        .splitCsv(header: true, sep: getSep(metadata_delim))
        .map { row ->
            if (!row.containsKey('sample_id') || !row.containsKey('group')) {
                error "Metadata missing required columns ('sample_id', 'group'). Found: ${row.keySet()}"
            }
            if (!row.sample_id || row.sample_id.trim() == "") {
                error "There are missing sample_id values in metadata"
            }
            if (!row.group || row.group.trim() == "") {
                error "There are missing group values in metadata"
            }
            return [ sample_id: row.sample_id, group: row.group, raw_row: row ]
        }

    // Parse & validate samplesheet
    ch_samplesheet = Channel
        .fromPath(samplesheet_path, checkIfExists: true)
        .splitCsv(header: true, sep: getSep(samplesheet_delim))
        .map { row ->
            if (!row.containsKey('input_path') || !row.containsKey('delimiter')) {
                error "Sample sheet missing required columns ('input_path', 'delimiter'). Found: ${row.keySet()}"
            }
            if (!['csv', 'tsv'].contains(row.delimiter)) {
                log.warn "WARNING: Unsupported delimiter '${row.delimiter}' for file ${row.input_path}"
            }
            return tuple(file(row.input_path, checkIfExists: true), row.delimiter)
        }

    emit:
    metadata    = ch_metadata.collect()
    samplesheet = ch_samplesheet.collect(flat: false)
}

workflow BUILD_SAMPLE_REGISTRY {
    take:
    ch_metadata     // Metadata channel
    ch_samplesheet  // Samplesheet channel

    main:
    // Flatten ch_samplesheet first to unpack the collected ArrayBag
    ch_sample_to_matrix = ch_samplesheet
        .flatten()
        .collate(2) // Groups the elements back into [file, delimiter] tuples
        .flatMap { file_obj, delim ->
            if (!['csv', 'tsv'].contains(delim.toLowerCase())) {
                return []
            }

            def headerLine = file(file_obj).withReader { reader -> reader.readLine() }
            if (!headerLine) return []

            def separator = getSep(delim)
            def columns = headerLine.split(separator)*.trim()

            // Filter out non-sample columns
            def sampleIds = columns.findAll { it && it != 'Mutation Types' }

            return sampleIds.collect { sample_id ->
                [ sample_id, file_obj.toString(), delim ]
            }
        }

    // Match extracted sample IDs against metadata
    ch_registry = ch_metadata
        .flatten()
        .map { meta ->
            [ meta.sample_id, meta ]
        }
        .join(ch_sample_to_matrix, by: 0)
        .map { sample_id, metadata, matrix_path, delim ->
            [
                sample_id  : sample_id,
                matrix_path: matrix_path,
                delim      : delim
            ]
        }

    emit:
    registry = ch_registry
}