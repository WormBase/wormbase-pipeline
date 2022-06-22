import os

def bash_check_if_success_or_loop(bcommand, times):
    if times < 2:
        times = 2
    if not bcommand.endswith("\n"):
        bcommand += "\n"
    bash_command = "for run in {1.." + str(times) + "}; do" + "\n" + \
                   "\t" + bcommand + \
                   "\t" + 'if [ "$?" = "0" ] ; then' + "\n" + \
                   "\t\t" + 'echo "completed normally";' + "\n" + \
                   "\t\t" + "break;" + "\n" + \
                   "\t" + "else" + "\n" + \
                   "\t\t" + 'echo "failed. Retrying...";' + "\n" + \
                   "\t\t" + "sleep 60;" + "\n" + \
                   "\t" + "fi" + "\n" + \
                   "done;\n"
    return (bash_command)


def md5_check(original_file, copied_file):
    bash_command = "\n\noriginal_md5=$(md5sum " + original_file + " | cut -d' ' -f1);\n" + \
                   "copied_md5=$(md5sum " + copied_file + " | cut -d' ' -f1);\n" + \
                   'if [ "$original_md5" != "$copied_md5" ];\n' + \
                   "\tthen echo 'File copy was unsuccessful. Exiting!'; exit 1;" + \
                   "fi;\n\n\n"
    return (bash_command)


def rsync(origin, destination):
    bash_command = "rsync -avz --checksum " + \
                   origin + " " + \
                   destination + ";\n"
    return (bash_command)


def delete_if_ok(todelete):
    bash_command = "if [ $? -eq 0 ]; then\n" + \
                   '\trm ' + todelete + '\n' + \
                   'else\n' + \
                   '\techo 1>&2 "command failed";\n' \
                   '\texit 1;\n' + \
                   'fi;\n\n'
    return bash_command


def do_if_exist(tfile, todo):
    bash_command = "if [ -f " + tfile + " ]; then\n" + \
                   todo + ';\n' + \
                   'else\n' + \
                   '\techo 1>&2 "File ' + tfile + ' does not exist.";\n' \
                                                  '\texit 1;\n' + \
                   'fi;\n\n'
    return bash_command


def expand_aliases():
    return ("shopt -s expand_aliases\n\n")
