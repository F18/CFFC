import sys, os, string
prefix = '/usr/local/cantera' 
pycmd = '/usr/bin/python2' # sys.argv[2]
localinst = 0

build_python = 0
build_matlab = 0

bindir = '/usr/local/cantera/bin'
libdir = '/usr/local/cantera/lib' 
hdrdir = '/usr/local/cantera/include/cantera'
demodir = '/usr/local/cantera/demos'
datadir = '/usr/local/cantera/data'
templdir = '/usr/local/cantera/templates'
ctdir = '/usr/local/cantera'
home = '/nfs/kris/d1/people/charestm'
isdarwin = '0'
numarray_home = ''

f = open(home+'/setup_cantera','w')
f.write('#!/bin/sh\n')
f.write('LD_LIBRARY_PATH='+libdir+':$LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH\n')
f.write('PATH='+bindir+':$PATH\nexport PATH\n')
f.write('PYTHON_CMD='+pycmd+'\nexport PYTHON_CMD\n')
if pycmd <> 'python':
    f.write('alias ctpython='+pycmd+'\n')

ctloc = '-'
warn = ''
warn2 = ''
pypath = ''
if localinst and build_python == 2:
    try:
        v = sys.version_info
        ctloc = prefix+'/lib/python'+`v[0]`+'.'+`v[1]`+'/site-packages'
        try:
            import Cantera
            ctpath = Cantera.__path__[0]
            if ctpath <> ctloc+'/Cantera':
                warn = """
 ######################################################################
    Warning: the Cantera Python package is already installed at
    """+ctpath+""".
    The newly-installed package at
    """+ctloc+"""/Cantera
    cannot be accessed until the existing one is removed.
 ######################################################################

"""
        except:
            pass

        pypath = ctloc
        if numarray_home <> '':
            pypath += ':'+numarray_home+'/lib/python'
            
        sys.path.append(ctloc)
        sys.path.append(numarray_home+'/lib/python')
        f.write('PYTHONPATH='+pypath+':$PYTHONPATH\nexport PYTHONPATH\n')
    except:
        print 'error'
f.close()


if localinst and build_python == 1:
    try:
        v = sys.version_info
        ctloc = prefix+'/lib/python'+`v[0]`+'.'+`v[1]`+'/site-packages'
        f.write('PYTHONPATH='+ctloc+':$PYTHONPATH\nexport PYTHONPATH\n')
    except:
        print 'error'
f.close()


if build_python == 2:
    # write the script to run MixMaster
    f = open(bindir+'/mixmaster','w')
    if isdarwin == '1':
        f.write('#!/bin/sh\n'+pycmd+"""w -c 'from MixMaster import MixMaster; MixMaster()'
        """)
    else:
        f.write('#!/bin/sh\n'+pycmd+""" -c 'from MixMaster import MixMaster; MixMaster()'
        """)
    f.close()


# write the script to copy files to build a new app
f = open(bindir+'/ctnew','w')
f.write("""#!/bin/sh
if test "x$1" = "x-f77"; then
cp """+templdir+"""/f77/*.* .
elif test "x$1" = "x-f90"; then
cp """+templdir+"""/f90/*.* .
else
cp """+templdir+"""/cxx/*.* .
fi
""")
f.close()


try:
    import Cantera
    ctpath = Cantera.__path__[0]
except:
    print "Cantera not found on sys.path = ",sys.path
    ctpath = "-"

if build_matlab:
    fm = open(ctdir+"/ctpath.m","w")
    fm.write("""path('"""+prefix+"""/matlab/toolbox/cantera/cantera',path)\n""")
    fm.write("""path('"""+prefix+"""/matlab/toolbox/cantera/cantera/1D',path)\n""")
    fm.close()

    fm = open(ctdir+"/cantera_demos.m","w")
    fm.write("""ctpath;\n""")
    fm.write("""cd demos/matlab;\n""")
    fm.write("""run_examples;\n""")
    fm.close()

print """

Cantera has been successfully installed.

File locations:

    applications      """+bindir+"""
    library files     """+libdir+"""
    C++ headers       """+hdrdir+"""
    demos             """+demodir+"""
    data files        """+datadir

if build_matlab:
    print """
    
    Matlab toolbox    """+prefix+"""/matlab/toolbox/cantera/cantera
    Matlab demos      """+prefix+"""/matlab/toolbox/cantera/cantera-demos
    Matlab tutorials  """+prefix+"""/matlab/toolbox/cantera/cantera-tutorials

    An m-file to set the correct matlab path for Cantera
    is at """+ctdir+"""/ctpath.m"""
    
if ctpath <> "-":
    print """
    Python package    """+ctpath
    if warn <> '':
        print warn
elif build_python == 2:
    print """
 ######################################################################
    Warning: the Cantera Python package is not installed. If you
    intentionally skipped it, ignore this message. Otherwise, type
    'make python' and/or 'make python-install' and look for error messages.
    Note that you must first install the 'numarray' package before installing
    the Cantera package.
 ######################################################################            
"""


print """
    
    setup script      """+home+"""/setup_cantera
    
    The setup script configures the environment for Cantera. It is
    recommended that you run this script by typing

      source """+home+"""/setup_cantera
    
    before using Cantera, or else
    include its contents in your shell login script.
    """

    
