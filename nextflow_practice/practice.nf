workflow {

    reads_ch = Channel.fromPath("test/*.fastq")

    reads_ch.view()

    // assembly_ch = Channel.fromPath("test/*.fasta")

    // делаем шаг один
    FASTQC(reads_ch)

    // делаем шаг два
    // MINIMAP2(input_ch)

    // делаем шаг три
    // input_ch = MINIMAP2.out.collect()
    // GET_COVERAGE(input_ch)

}


process FASTQC {
    conda 'envs/fastqc_env.yaml'
    container 'biocontainers/fastqc:v0.11.9_cv7'


    input:
    path(input_file)
    
    output:
    path("fastqc_output/*.html")
    
    script:
    """
    mkdir fastqc_output && fastqc ${input_file} -o fastqc_output/
    """
}

process MINIMAP2 {
    input:
    path(input_file)
    
    output:
    path("result.txt")
    
    script:
    """
    process.sh ${input_file} > result.txt
    """
}

process GET_COVERAGE {
    input:
    path(input_file)
    
    output:
    path("result.txt")
    
    script:
    """
    process.sh ${input_file} > result.txt
    """
}