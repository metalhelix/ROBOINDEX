process faidx_Make {

    label 'big_mem'

    input: 
    val meta

    output:
    path "*.fai"

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:
    def fna = file(meta.fna_path)
    """
    ml samtools

    samtools faidx ${fna} 

    mv ${meta.fna_path}.fai ./
    """
}

process bt2Index_Make {

    label 'big_mem'

    input:
    val meta

    output:
    path "bowtie2/*" //The naming might be different if the fasta file does not end with .fa

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:
    """
    ml bowtie2

    mkdir -p bowtie2

    bowtie2-build ${meta.fna_path} bowtie2/${meta.id}
    """
}

process bwaIndex_Make {
    
    label 'big_mem'

    input:
    val meta

    output:
    path "BWA/*"

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:

    def fna = file(meta.fna_path)

    """
    ml bwa

    mkdir -p BWA

    bwa index ${fna}

    mv `ls ${fna}.*|grep -v ".fai"|grep -v ".dict"` BWA
    """
}

process cellranger_singlecell_Index_Make {

    memory { 40.GB * task.attempt }
    time { 10.hour * task.attempt }
    cpus { 10 }
    
    errorStrategy 'retry'
    maxRetries 3

    input:
    val meta

    output:
    path "10x/*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/", mode: 'copy'

    script:
    """
    ml cellranger/8.0.1 
    
    cellranger_ver=cellranger_8_0_1_proteinCoding 

    mkdir -p 10x 

    cellranger mkgtf ${meta.gtf_path} filtered.gtf --attribute=gene_biotype:protein_coding

    cellranger mkref --genome=\$cellranger_ver  --fasta=${meta.fna_path} --genes=filtered.gtf --memgb=300 

    ln -s \$cellranger_ver cellranger

    mv \$cellranger_ver cellranger filtered.gtf 10x
    """
}

process igv_genome_Make {
    
    label 'lil_mem'
    
    input: 
    val meta 
    val riboList_Make_out //tells the pipeline to wait until that process is finished

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/extras", mode: 'copy'
    
    output:
    path "${meta.id}.${meta.annotation_version}.igv_genome.json"

    script:
    """
    igv_genome_make_v2.py -g ${meta.index_Folder}/annotation/${meta.annotation_version}/gtfs/${meta.id}.${meta.annotation_version}.gtf -n ${meta.index_Folder}/${file(meta.fna_path).getName()} -f ${meta.index_Folder}/${file(meta.fna_path).getName()}.fai -i ${meta.id} -a ${meta.name} -v ${meta.annotation_version}
    """
}

process hisat2_index_build {

    label 'big_mem'

    input:
    val meta

    output:
    path "*.ht2"
    path "hisat2_version.txt"
    
    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/HISAT2", mode: 'copy'
    
    script:
    """
    ml hisat2

    hisat2-build ${meta.fna_path} ${meta.id}

    hisat2-build --version > hisat2_version.txt
    """
    }
