import os

def rsync_to_external_location_command(directory_from, remote_path, method="embassy"):
    if method == "embassy":
        return(rsync_dir_to_external_location_embassy_bucket_command(directory_from, remote_path))
    else:
        exit_with_error("You can only copy files to embassy right now.")

def rsync_dir_to_external_location_embassy_bucket_command(directory_from, embassy_path):
    bash_command = "module load embassy\n\n" + \
                   "eaws s3 sync " + directory_from + " " + \
                   embassy_path
    return bash_command