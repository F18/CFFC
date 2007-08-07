import sys
try:
    from distutils.core import setup, Extension
except:
    print 'could not import distutils. Will try anyway...'
    
libs = []
platform = sys.platform

flibs = ''
#linkargs = ''

numarray_incl = ""

incdirs=["../../build/include", "src", "../clib/src"]

if numarray_incl <> '':
    incdirs.append(numarray_incl)


bllibstr = "-lctlapack -lctblas"
bllibs = bllibstr.replace('-l',' ')
bllist = bllibs.split()

cvlibstr = "-lcvode"
cvlibs = cvlibstr.replace('-l',' ')
cvlist = cvlibs.split()

extra_link = ""
linkargs = extra_link.split()    

bldirstr = " -L/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu"
bldirs = bldirstr.replace('-L',' ')
dirlist = bldirs.split()
libdir = ['/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu']
for d in dirlist:
    libdir.append(d)

endlibstr1 = " -lctf2c -lm"
endlib1 = endlibstr1.replace('-l', ' ')
endlib = endlib1.split()

if platform == "win32":
    libs = ["clib", "zeroD","oneD","transport",
            "cantera","recipes"] + bllist + cvlist + ["ctmath", "tpx", "converters"]
else:
    
    libs = ["clib", "zeroD","oneD","transport",
            "cantera", "converters"] + bllist + cvlist + ["ctmath", "tpx"]
                                          
    if 1 == 1:
         libs.append("ctf2c")
    for d in endlib:
          libs.append(d) 


# values:
#  0   do nothing
#  1   install only ctml_writer.py
#  2   install full package
#  3   try to install full, but install ctml_writer if full package
#      install fails
buildPython = 0
if buildPython >= 2:

    try:
        setup(name="Cantera",
              version="1.7.0",
              description="The Cantera Python Interface",
              long_description="""
              """,
              author="Prof. D. G. Goodwin, Caltech",
              author_email="dgoodwin@caltech.edu",
              url="http://www.cantera.org",
              package_dir = {'MixMaster':'../../apps/MixMaster'},
              packages = ["","Cantera","Cantera.OneD",
                          "MixMaster","MixMaster.Units"],
              ext_modules=[ Extension("Cantera._cantera",
                                      ["src/pycantera.cpp"],
                                      include_dirs=incdirs,
                                      library_dirs = libdir,
                                      libraries = libs,
                                      extra_link_args = linkargs
                                      )
                            ],
              )
    except:
        buildPython = 1
    
        
if buildPython == 1:
    try:
        setup(name="Cantera CTI File Processor",
              version="1.7.0",
              description="Converts .cti files to CTML",
              long_description="""
              """,
              author="Prof. D. G. Goodwin, Caltech",
              author_email="dgoodwin@caltech.edu",
              url="http://www.cantera.org",
              py_modules = ["ctml_writer"],
              )
    except:
        raise 'Error encountered while building or installing the Cantera CTI file preprocessor!'
