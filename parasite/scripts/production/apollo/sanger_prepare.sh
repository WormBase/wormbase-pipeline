function sanger_prepare {
    ssh -f -N -C -X ssh.sanger.ac.uk
    ssh -f -N -C -X sangerfarm
}

function sanger_connect {
    sanger_prepare;
    ssh sangerngs;
}



