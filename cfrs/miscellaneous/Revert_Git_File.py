usage = """%prog [options]

DESCRIPTION:
     Reverts a FILE to the last committed revision.
"""

import re, commands, sys, optparse, os
from optparse import OptionParser

# ============================= INPUT PARSING ============================
class OptionParserNoError(OptionParser):
    def error(self, msg):
        return 

def ParseCommandLine():
    ####################### Parse the command line ######################
    parser = OptionParser(usage=usage)
    
    parser.add_option('-f', '--file', action="store",
                      dest='filename', metavar="FILE", type='string',
                      help='specify the error norms input file')
   
    (options, args) = parser.parse_args()

    if (options.filename == None):
        parser.error("An input file must be provided!!!\n" + "Use --help to see the options.")
    elif os.path.isfile(options.filename) == False:
        parser.error(options.filename + " is not a valid input file")

    # Return
    return options, args


# ==========================================
# ******************   MAIN  ***************
# ==========================================

#### Parse the command line #######
(options, args) = ParseCommandLine()

# Target file --> the file that is finaly modified
Target = options.filename

################# Test for existance of 'git'
cmd = "which git"
cmdOutput = commands.getoutput(cmd)       # get the output of this command

GitExist = cmdOutput.find('no git in')

if (GitExist != -1) :                       # Git is not installed on the system
    sys.exit(0)
else:
    cmd = 'git checkout ' + options.filename
    cmdOutput = commands.getoutput(cmd)       # get the output of this command
### DONE ####
