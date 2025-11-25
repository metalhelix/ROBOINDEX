process SampleSheet_SourceEnsembl {
    
    memory { 10.GB * task.attempt }
    time { 5.hour * task.attempt }
    cpus { 4 }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path samplesheet

    output:
    path "samplesheet_processed.csv"

    script:
    """
    ensembl_fetch.py -s ${samplesheet} -d ${params.out}
    """
}