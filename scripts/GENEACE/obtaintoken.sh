#!/bin/bash

client_id=$1
client_secret=$2

USAGE="obtaintoken.sh \$CLIENT_ID \$CLIENT_SECRET"

if [[ -z "$client_id" || -z "$client_secret" ]]
then
    echo "$USAGE"
    exit 1
fi

test $(uname -o) = "GNU/Linux"
if [ $? -eq 0 ]; then
    OPEN=xdg-open
else
    OPEN=open
fi

redirect_uri="urn:ietf:wg:oauth:2.0:oob"
scope="email%20profile%20openid"
$OPEN "https://accounts.google.com/o/oauth2/auth?client_id=${client_id}&redirect_uri=$redirect_uri&scope=$scope&response_type=code" \
      2> /dev/null 3> /dev/null

read -sp "Paste Google one-time-code:" code
RESP=$(curl --silent \
	    -d "code=$code" \
	    -d "client_id=$client_id" \
	    -d "client_secret=$client_secret" \
	    -d "redirect_uri=$redirect_uri" \
	    -d "grant_type=authorization_code" \
	    https://accounts.google.com/o/oauth2/token)

echo $RESP | jq .id_token | tr -d '"'
