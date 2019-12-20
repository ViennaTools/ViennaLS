#!/bin/bash

#root directory of project and config file for doxygen
parentname="ViennaLS"
templatefile="doxygen_template.txt"
configfile="doxygen_config.txt"

#filename for html shortcut to Documentation outside of html folder
docmainfile="../index.html"


#ViennaLS directory is automatically found by searching current working directory
viennaLSdir=$(pwd | sed 's/\(.*ViennaLS\).*/\1/' | sed -e 's/[\/&]/\\&/g')

echo "Using directory: ${viennaLSdir}"
#list of all commands and arguments that need to be changed for different computers
commands[0]="USE_MDFILE_AS_MAINPAGE"
commands[1]="INPUT"
commands[2]="EXCLUDE"
commands[3]="PROJECT_LOGO"
commands[4]="STRIP_FROM_PATH"
commands[5]="EXAMPLE_PATH "

argument[0]="$viennaLSdir\/README.md"
argument[1]="$viennaLSdir" 
argument[2]="$viennaLSdir\/cmake $viennaLSdir\/Tests $viennaLSdir\/Wrapping"
argument[3]="$viennaLSdir\/docs\/doxygen\/logo.png"
argument[4]="$viennaLSdir"
argument[5]="$viennaLSdir\/Examples"

#create doxygen config file
cp $templatefile $configfile

#change the doxygen config file to the new parameters
for i in `seq 0 $(( ${#commands[*]} - 1 ))`; do
	command='sed -i "s/\(.*'${commands[$i]}' *=\).*/\1 '${argument[$i]}'/" '$configfile
	eval $command
done

#run doxygen to make the page
doxygen $configfile

#if shortcut outside the html folder does not exist, create it
if [ -e $docmainfile ]
then
	echo "$docmainfile already exists: Not creating."
else
echo '<!DOCTYPE html>
<html lang="en-US">
    <head>
        <meta charset="UTF-8">
        <meta http-equiv="refresh" content="1; url=doxygen/html/index.html">
        <script type="text/javascript">
            window.location.href = "doxygen/html/index.html"
        </script>
        <title>Page Redirection</title>
    </head>
    <body>
        If you are not redirected automatically, follow this <a href="doxygen/html/index.html">link to the manual</a>.
    </body>
</html>' > $docmainfile
echo "Created shortcut: $docmainfile"
fi

