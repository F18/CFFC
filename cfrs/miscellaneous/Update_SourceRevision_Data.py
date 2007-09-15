import re
import commands

# Target file --> the file that is finaly modified
Target = "Common/SourceRevisionData.cc"

# Dictionaries
# RevisionData --> associated with GIT output
RevisionData = {
    'CompilationTime' : ' ',
    'Last Changed Author' : ' ',
    'Last Changed Commit' : ' ',
    'Last Changed Date': ' '
}

# CodeRevisionData --> associated with the '*.cc' file. That's the final modified file
CodeRevisionData = {
    'SourceCode::CompilationTime' : ' ',
    'SourceCode::LastChanged_Author' : ' ',
    'SourceCode::LastCommitted_Hash' : ' ',
    'SourceCode::LastChanged_Date' : ' ',
}

# Translator --> describes the conversion between the keys of RevisionData and CodeRevisionData
Translator = {
    'CompilationTime'     : 'SourceCode::CompilationTime',
    'Last Changed Author' : 'SourceCode::LastChanged_Author',
    'Last Changed Date'   : 'SourceCode::LastChanged_Date',
    'Last Changed Commit'    : 'SourceCode::LastCommitted_Hash',
}

################# Get compilation time with 'date' ###############
cmd="date '+%d %b %Y, %R'"
RevisionData['CompilationTime'] = commands.getoutput(cmd)

################# Get GIT info ###################
cmd="git log --pretty=format:\"Last Changed Commit= %H %nLast Changed Author= %an %nLast Changed Date= %cd %n\" -1"
gitinfo = commands.getoutput(cmd)       # get the output of this command

####### Update the values for RevisionData
list = gitinfo.split("\n")

for entry in list:
    new_entry = entry.strip().split("=")

    for key in RevisionData.keys():
        if new_entry[0] == key:
            if key == 'Last Changed Date':
                date = new_entry[1].strip().split(' ')
                # format the 'Last Changed Date'
                RevisionData[key] = date[2] + " " + date[1] + " " + date[4] + ", " + date[3][0:5]
            else :
                RevisionData[key] = new_entry[1].strip()

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
