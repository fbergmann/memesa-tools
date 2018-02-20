"""
Convert SBML model description to the H-presentation (.ine file)
===============================================================

This script accepts command line arguments (-m *model_name*  -b *inf_bound*  -l *sbml_level*) but you can also specify these arguments in this script

User input:
- *model_name* (string)
- *inf_bound* (integer/float)
- *sbml_level* (integer)

This scripts accepts both SBML L2 and L3 input files, but L3 are preferred.


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

(C) B.G. Olivier, Timo R. Maarleveld, M.T. Wortel, F.J. Bruggeman, and B. Teusink

Written by B.G. Olivier Amsterdam, TR Maarleveld, The Netherlands
E-mail: b.g.olivier@vu.nl
Last Change: Dec 18, 2017
"""

from __future__ import division, print_function, absolute_import

### Set userdata here ####

model_name = 'toy_model'
inf_bound = 1000
sbml_level = 2

### End userdata ####

import os, time, numpy, sys,getopt
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

try:
    myopts, args = getopt.getopt(sys.argv[1:],"m:b:l:")
except getopt.GetoptError as e:
    print(str(e))
    print("Usage: %s -m model_name -b inf_bound -l sbml_level" % sys.argv[0])
    sys.exit(2)

###################################
# opt == option
# arg == argument passed to the opt
###################################
for opt, arg in myopts:
    if opt == '-m':
        model_name = arg
    elif opt == '-l':
        sbml_level = int(arg)
    elif opt == '-b':
        inf_bound = float(arg)
    else:
        print("Usage: %s -m model_name -b inf_bound -l sbml_level" % sys.argv[0]) 

if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')
    
    
# bgoli 2017-06-21
# CBMPy compatability changes
import cbmpy
__DEBUG__ = cbmpy.CBConfig.__CBCONFIG__['DEBUG']
__version__ = cbmpy.CBConfig.current_version()
import cbmpy.fluxmodules.fluxmodules as fmod
from cbmpy import CBRead, CBWrite, CBTools,CBMultiCore
from cbmpy import CBSolver as slv

### Active script ###

model_file = model_name + '.xml'
data_dir = os.path.join(cDir, 'data', model_name)
sbml_dir = os.path.join(data_dir,'models_subnetwork','sbml')
output_dir = os.path.join(data_dir,'cope_fba','subnetworks','vertex')
index_dir = os.path.join(data_dir,'models_subnetwork','h-format')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for f in os.listdir(sbml_dir):
    if 'correct' in f:
        submod = None
        idmap = {}
        F = open(os.path.join(index_dir, f.replace('.correct.xml', '.noinf_r.columns.txt')), 'r')
        for l in F:
            L = l.split(',')
            print(L)
            if len(L) == 2:
                idmap[L[1].strip()] = int(L[0])+1
        F.close()
        print(idmap)
        if sbml_level == 3:
            print('Trying level 3')
            submod = CBRead.readSBML3FBC(f, work_dir=sbml_dir)
        else:
            submod = CBRead.readSBML2FBA(f, work_dir=sbml_dir)
        if submod is not None:
            cbmpy.doFBA(submod)
            a,b = cbmpy.doFVA(submod)
            output = {}
            for r in submod.reactions:
                span = abs(r.fva_max - r.fva_min)
                rtype = 'VARIABLE'
                if span < 1e-10:
                    rtype = 'FIXED'
		if numpy.isinf(float(r.fva_max)): #EvP changed because of numpy.inf is not rational
		    fva_max = float(999999)
		else:
		    fva_max = r.fva_max
                S = 'x{} : {} : {} -- {}'.format(idmap[r.getId()], r.fva_min, fva_max, rtype)
                print(S)
                output[idmap[r.getId()]] = S 
            if len(output) > 0:
                kdx = list(output.keys())
                kdx.sort()
                outstr = ""
                for i in kdx:
                    outstr += '{}\n'.format(output[i])
                print(outstr)
                F = open(os.path.join(output_dir, f.replace('.correct.xml', '.noinf_r.ine.opt.fva')), 'w')
                F.write(outstr)
                F.close()
                    
                
            
                


