include {GENOME_INDEX_BUILD; SECUNDO_INDEX_BUILD; BACTERIA_INDEX_BUILD} from "./workflows/Genome_Index.nf" // The main workflow
include {SAMPLESHEET_CHECK} from "./workflows/sampleSheet_check.nf" // To parse the sample sheet into the pipelin e
include {SampleSheet_SourceEnsembl} from "./modules/fetch-ensembl_process.nf"

workflow ROBO_FULL {
    take:
    samplesheet

    main:
    SampleSheet_SourceEnsembl(samplesheet)
    SAMPLESHEET_CHECK( SampleSheet_SourceEnsembl.out)
    meta = SAMPLESHEET_CHECK.out.data_meta
    GENOME_INDEX_BUILD(meta)
}

workflow ROBO_SECUNDO {
    take:
    samplesheet

    main:
    SampleSheet_SourceEnsembl(samplesheet)
    SAMPLESHEET_CHECK( SampleSheet_SourceEnsembl.out) 
    meta = SAMPLESHEET_CHECK.out.data_meta
    SECUNDO_INDEX_BUILD(meta)
}

workflow ROBO_BACTERIA {
    take:
    samplesheet

    main:
    SampleSheet_SourceEnsembl(samplesheet)
    SAMPLESHEET_CHECK( SampleSheet_SourceEnsembl.out) 
    meta = SAMPLESHEET_CHECK.out.data_meta
    BACTERIA_INDEX_BUILD(meta)
}


if (!file(params.samplesheet).exists()) {
    exit 1, "ERROR: Please check input samplesheet -> ${params.samplesheet}"
} else if (!params.samplesheet) { 
    exit 1, "ERROR: Please provide a samplesheet "
} else {
    workflow {
        if (!params.secundo) {
            if (!params.bacteria) {
                ROBO_FULL(params.samplesheet)
                } else {
                ROBO_BACTERIA(params.samplesheet)
                }
        } else {
            ROBO_SECUNDO(params.samplesheet)
        }
    }
}
