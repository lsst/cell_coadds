#!/usr/bin/bash

# Source: https://auscunningham.medium.com/enforcing-git-commit-message-style-b86a45380b0f
# Source: https://git-scm.com/book/en/v2/Customizing-Git-An-Example-Git-Enforced-Policy

commit_message_check(){
    #  Get all the commit hashes in this branch
    git rev-list HEAD...origin/master > .shalist.tmp.txt
    for i in `cat .shalist.tmp.txt`
    do
        gitmessage=`git cat-file commit "$i" | sed '1,/^$/d'`
        header=`echo "$gitmessage" | head -n 1`

        # Check that the first line begins with an upper case letter
        if [[ ! $header =~ ^[A-Z] ]]
        then
            echo "[ERROR] Commit message does not begin with an upper case letter: $header ($i)"
        fi

        # Check that the first line is within 50 characters
        if [ ${#header} -gt 50 ]
        then
            echo "[ERROR] Commit message header too long: $header ($i)"
        fi

        # Check that the first line does not end with a period
        if [[ $header =~ \.$ ]]
        then
            echo "[ERROR] Commit message header ends with a period: $header ($i)"
        fi

        # Chech that there are no trailing whitespaces in the header
        if [[ $header =~ \s$ ]]
        then
            echo "[ERROR] Commit message header ends with a whitespace: $header ($i)"
        fi

        # Warn if there is and or & in the commit message header
        if [[ $header =~ \& ]]
        then
            echo "[WARNING] Commit message header contains an ampersand: $header ($i)"
        fi

        # Check that the second line is empty (or non-existent)
        if [ "$(echo "$gitmessage" | head -n 2 | tail -n 1)" != "" ]
        then
            echo "[ERROR] Commit message has a second line: $header ($i)"
            exit 1
        fi

        # Fetch the body of the commit message
        body=`echo $gitmessage | sed '1,/^$/d' | sed '1,/^$/d'`
        if [ ! -z "$body" ]
        then
            echo "Body is present"
            exit 0
        fi


    done
    rm .shalist.tmp.txt
}

# Call the function
commit_message_check