def trim_paired_fqs(fastp, input_fqs, output_fqs,
                    fastp_html_report, fastp_json_report,
                    threads):
    ispaired = 1
    if isinstance(input_fqs, list) and len(input_fqs) == 1:
        ispaired = 0
    elif isinstance(input_fqs, str):
        ispaired = 0
        input_fqs = [input_fqs]
        output_fqs = [output_fqs]
    bash_command = fastp + " -i " + input_fqs[0] + " " + \
                   ("-I " + input_fqs[1] if ispaired==1 else "") + " " + \
                   "-o " + output_fqs[0] + " " + \
                   ("-O " + output_fqs[1] if ispaired==1 else "") + " " + \
                   "--html " + fastp_html_report + " " + \
                   "--json " + fastp_json_report + " " + \
                   "-w " + str(threads) + ";\n\n"
    return(bash_command)



def curl_fastqs(fastqs_dict):
    """fastqs_dict should be a dictionary with the remote paths as keys and the output paths as values"""
    bash_command = "\n\n".join(
        ["curl -H 'Connection: keep-alive' --keepalive-time 2 -o " + fastqs_dict[fastq] + " " + fastq + ";" for fastq in fastqs_dict]) + "\n\n"
    return (bash_command)