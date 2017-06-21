"""
Find Duplicate Reactions
========================

Please remove all modelling duplicates, artefacts and  isozymes.

Biological relevant ones should stay in the model (e.g. atp maintenance vs nucleotide biosynthesis)

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

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, F.J. Bruggeman, and B. Teusink

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last Change: February 25, 2015
"""
from __future__ import division, print_function, absolute_import

# Set userdata here #

model_name = 'toy_model'
inf_bound = 1000
sbml_level = 3

# End userdata #

import sys,getopt

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
        print("Usage: %s -i input -o output" % sys.argv[0]) 
        
### Active script ###

if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')

import os, time, numpy
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
model_file = model_name + '.xml'
data_dir = os.path.join(cDir, 'data', model_name)
model_dir = os.path.join(data_dir,'models','sbml')
work_dir = os.path.join(data_dir,'models','h-format')

# bgoli 2017-06-21
# CBMPy compatability changes
import cbmpy
__DEBUG__ = cbmpy.CBConfig.__CBCONFIG__['DEBUG']
__version__ = cbmpy.CBConfig.current_version()

from cbmpy import CBRead, CBWrite, CBTools
from cbmpy import CBSolver as slv

work_dir = os.path.join(cDir, work_dir)
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

if sbml_level == 3:
    print('Trying level 3')
    cmod = CBRead.readSBML3FBC(model_file, work_dir=model_dir)
else:
    cmod = CBRead.readSBML2FBA(model_file, work_dir=model_dir)
cmod.id = model_name

slv.analyzeModel(cmod, lpFname=os.path.join(work_dir, 'raw_({0:s})'.format(cmod.id)))

# scan for duplicates
print('\nScan for duplicates')
dup_C = CBTools.scanForReactionDuplicates(cmod, ignore_coefficients=False)
print('\nReaction pairs with matching reagents and coefficients: {0:d}'.format( len(dup_C) ))
for d in dup_C:
    print('{0:s} == {1:s} ({2:s})'.format(d[0],d[1],d[2]))

# write duplicate reactions to file
F = open(os.path.join(model_dir,'{0:s}.duplicates.csv'.format(model_file)),'w')
F.write('Model file: {0:s},,,,\n'.format( model_file) )
F.write('Duplicate reactions pairs: {0:d},,,,\n'.format( len(dup_C) ) )
for d in dup_C:
    F.write('\"{0:s}\",\"{1:s}\",\"{2:s}\",\"{3:s}\",\"{4:s}\"\n'.format(d[0],d[1],d[2],d[3],d[4]))
F.write('\n')
F.close()
