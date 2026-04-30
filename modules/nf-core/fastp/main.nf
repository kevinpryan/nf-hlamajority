process FASTP {
    tag "$meta.sample"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0':
        'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)
    path  adapter_fasta
    val   save_trimmed_fail
    val   save_merged

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    
    // --- DYNAMIC LOGIC START ---
    // Check if Single End (either via meta flag or list size)
    def single_end = meta.single_end || !(reads instanceof List)
    // 1. Prepare Inputs
    def fastp_in  = "--in1 ${prefix}_1.fastq.gz"
    def symlink_cmds = "[ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz"
    
    if (!single_end) {
        fastp_in += " --in2 ${prefix}_2.fastq.gz"
        symlink_cmds += "\n[ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz"
    }

    // 2. Prepare Outputs
    def fastp_out = "--out1 ${prefix}_1.fastp.fastq.gz"
    if (!single_end) {
        fastp_out += " --out2 ${prefix}_2.fastp.fastq.gz"
    } else {
        // For SE, we usually just want the output named appropriately (e.g. without _1 if you prefer, but keeping _1 makes regex easier)
        // Leaving it as _1.fastp.fastq.gz is fine and keeps wildcard *.fastp.fastq.gz working
    }

    // 3. Prepare "Fail" logs
    def fail_fastq = ""
    if (save_trimmed_fail) {
        if (!single_end) {
            fail_fastq = "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz"
        } else {
            fail_fastq = "--failed_out ${prefix}.fail.fastq.gz"
        }
    }

    // 4. Prepare Merging (Only valid for PE)
    def merge_fastq = (!single_end && save_merged) ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
    
    // 5. PE specific flags
    def pe_flags = (!single_end) ? "--detect_adapter_for_pe" : ""
    // --- DYNAMIC LOGIC END ---

    """
    $symlink_cmds

    fastp \\
        $fastp_in \\
        $fastp_out \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        $adapter_list \\
        $fail_fastq \\
        $merge_fastq \\
        --thread $task.cpus \\
        $pe_flags \\
        $args \\
        2> >(tee ${prefix}.fastp.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
