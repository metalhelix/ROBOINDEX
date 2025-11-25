
process STAR_Bacteria {
    memory { 30.GB * task.attempt }
    time { 10.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val meta 
    val cpus

    output:
    path "STAR*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}", mode: 'copy'

    script:
    """
    ml star 

    STAR --version > STAR.version
    
    STAR --runThreadN ${cpus} \
        --runMode genomeGenerate --genomeDir STAR \
        --genomeFastaFiles ${meta.fna_path} \
        --limitGenomeGenerateRAM=195729992629 \
        --genomeSAindexNbases 11

    """
}


