import json
import paramiko
import shutil
import os
import tempfile
from ProductionUtils import check_remote_file_exists,rename_remote_file

def write_to_remote_apollo_trackList_json(json_data, species, remote_address, write=True):

    # Save data to a temporary file
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        json.dump(json_data, temp_file, indent=4)
        temp_file_path = temp_file.name
    # Copy the temporary file to the remote location
    username, hostname, filepath = remote_address.split('@')[0], remote_address.split('@')[1].split(':')[0], remote_address.split(':')[1]

    destination_path = f"{username}@{hostname}:{filepath}"

    # Check if the remote file already exists
    if check_remote_file_exists(remote_address):
        # Rename the existing file with a .bak extension
        bak_filepath = f"{filepath}.bak"
        bak_destination_path = f"{username}@{hostname}:{bak_filepath}"
        if write:
            rename_remote_file(remote_address, bak_destination_path)


    # Copy the temporary file to the remote location
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(hostname=hostname, username=username)
    # Create an SFTP client
    sftp = ssh.open_sftp()
    try:
        # Copy the temporary file to the remote location
        if write:
            sftp.put(temp_file_path, filepath)
            print("File copied to remote location:", remote_address + filepath)
    finally:
        # Close the SFTP session and SSH connection
        sftp.close()
        ssh.close()

    # Optionally, return the destination path
    return destination_path

