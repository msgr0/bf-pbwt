params.exe = "/Users/msgro/Progetti/bf-pbwt/2bfpbwt"
params.output = ""
params.input = ""

process BFPBWT{
    storeDir "$params.output"

    input:
    tuple path(input), val(type)

    output:
    tuple path(time), path(out)

    script:
    time = "${input}.${type}.time.txt"
    out = "${input}.${type}.out.txt"
    """
    #!/usr/bin/env bash

    ${params.time_exe} -o ${time} ${params.exe} $type $input 2> ${out}
    """
}

workflow {
    // exe = Channel.fromPath("$params.exe")
    runtype_ch = Channel.of("bli", "bar", "blis", "bars", "bpr") | view
    input_ch = Channel.fromPath("$params.input/*.bm") | view


    input_ch = input_ch.combine(runtype_ch)

    BFPBWT(input_ch) | view



}
