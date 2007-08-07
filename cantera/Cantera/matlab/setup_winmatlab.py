
import sys

bindir = '/usr/local/cantera/bin'
libdir = '/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu'
incdir = '/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/include'
dflibdir = ''

bllibstr = "-lctlapack -lctblas"
bllibs = bllibstr.replace('-l',' ')
bllist = bllibs.split()

bldir = "/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu"

libs = ['clib', 'oneD', 'zeroD', 'transport', 'cantera', 'recipes',
        'cvode', 'ctmath', 'tpx']

f = open('setup.m','w')
f.write('cd cantera\nbuild_cantera\nexit\n')
f.close()

fb = open('cantera/build_cantera.m','w')
fb.write("""
disp('building Cantera..');
mex -I"""+incdir+""" private/ctmethods.cpp private/ctfunctions.cpp ...
    private/xmlmethods.cpp private/phasemethods.cpp  ...
    private/thermomethods.cpp private/kineticsmethods.cpp ...
    private/transportmethods.cpp private/reactormethods.cpp ...
    private/reactornetmethods.cpp ...
    private/wallmethods.cpp private/flowdevicemethods.cpp ...
    private/funcmethods.cpp ...
    private/onedimmethods.cpp private/surfmethods.cpp private/write.cpp ...
""")
s = ''
for lib in libs:
    s += '    '+libdir+'/'+lib+'.lib ...\n'
fb.write(s)
if bllist:
    s = ''
    for lib in bllist:
        s += '    '+bldir+'/'+lib+'.lib ...\n'
    fb.write(s)
fb.write('    "'+dflibdir+'/dformd.lib" ...\n')
fb.write('    "'+dflibdir+'/dfconsol.lib" ...\n')
fb.write('    "'+dflibdir+'/dfport.lib" \n')

fb.close()

fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

