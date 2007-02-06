#!/usr/bin/python

import re
import commands

# Target file --> the file that is finaly modified
Target = "Common/SourceRevisionData.cc"

# Dictionaries
# RevisionData --> associated with SVN output
RevisionData = {
    'CompilationTime' : ' ',
    'GlobalRevision': ' ',
    'Last Changed Author' : ' ',
    'Last Changed Rev' : ' ',
    'Last Changed Date': ' '
}

# CodeRevisionData --> associated with the '*.cc' file. That's the final modified file
CodeRevisionData = {
    'SourceCode::CompilationTime' : ' ',
    'SourceCode::LastChanged_Author' : ' ',
    'SourceCode::LastChanged_Revision' : ' ',
    'SourceCode::LastChanged_Date' : ' ',
    'SourceCode::Revision' : ' ',
}

# Translator --> describes the conversion between the keys of RevisionData and CodeRevisionData
Translator = {
    'CompilationTime'     : 'SourceCode::CompilationTime',
    'GlobalRevision'      : 'SourceCode::Revision',
    'Last Changed Author' : 'SourceCode::LastChanged_Author',
    'Last Changed Date'   : 'SourceCode::LastChanged_Date',
    'Last Changed Rev'    : 'SourceCode::LastChanged_Revision',
}

################# Get global revision with 'svnversion' ###############
cmd="svnversion"
RevisionData['GlobalRevision'] = commands.getoutput(cmd)

################# Get compilation time with 'date' ###############
cmd="date +%c"
RevisionData['CompilationTime'] = commands.getoutput(cmd)

################# Get SVN info ###################
cmd="svn info"
svninfo = commands.getoutput(cmd)       # get the output of this command

####### Update the values for RevisionData
list = svninfo.split("\n")
for entry in list:
    new = entry.strip().split(":")
    
    # for 'Last Changed Date' modify the second entry
    if new[0] == 'Last Changed Date':
        # match for the date between parenthesis
        match = re.search('\(.*\)',new[-1])
        date = match.group()[1:-1]
        # assign the date to the second entry of new
        new[1] = date
    
    for key in RevisionData.keys():
        if new[0] == key:
            RevisionData[key] = new[1].strip()

####### Update the new values for CodeRevisionData
for key in Translator:
    CodeRevisionData[Translator[key]] = RevisionData[key]

### Read Target file content
input = file(Target, "r")   # Open Target for reading
content = input.read()
input.close()

####### Update content of Target file
for key in CodeRevisionData.keys():
    input_pattern = key + "=\".*\""
    output_pattern = key + "=\"" + CodeRevisionData[key] + "\""
    content = re.sub(input_pattern,output_pattern,content)

### Write new content to Target file
output = file(Target, "w")   # Open Target for writing
output.write(content)
output.close()

### DONE ####
