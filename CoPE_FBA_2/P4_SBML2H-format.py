"""
Convert SBML model description to the H-presentation (.ine file)
===============================================================

This script accepts command line arguments (-m *model_name*  -b *inf_bound*  -l *sbml_level*) but you can also specify these arguments in this script

User input:
- *model_name* (string)
- *inf_bound* (integer/float)
- *sbml_level* (integer)

This scripts accepts both SBML L2 and L3 input files, but L3 are preferred.

Script adapted from the MEMESA-TOOLS project.

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

Written by TR Maarleveld, B.G. Olivier Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 27, 2015
"""

from __future__ import division, print_function, absolute_import

### Set userdata here ####

model_name = 'toy_model'
inf_bound = 1000
sbml_level = 3

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
     
from cbmpy.CBVersion import __DEBUG__, __version__
from cbmpy import CBRead, CBWrite, CBTools
from cbmpy import CBSolver as slv

### Active script ###

model_file = model_name + '.xml'
data_dir = os.path.join(cDir, 'data', model_name)
model_dir = os.path.join(data_dir,'models','sbml')
work_dir2 = os.path.join(data_dir,'models','h-format')

if not os.path.exists(work_dir2):
    os.makedirs(work_dir2)

if sbml_level == 3:
    print('Trying level 3')
    cmod = CBRead.readSBML3FBC(model_file, work_dir=model_dir)
else:
    cmod = CBRead.readSBML2FBA(model_file, work_dir=model_dir)

cmod.id = model_name
print('\nAttempting to delete bounds for biomass reaction,', cmod.getActiveObjective().getFluxObjectiveReactions()[0])
cmod.deleteBoundsForReactionId(cmod.getActiveObjective().getFluxObjectiveReactions()[0])

mLP = slv.analyzeModel(cmod, return_lp_obj=True)
CBWrite.printFBASolution(cmod)
tmp_mid = cmod.id+'_cplex'

CBTools.countedPause(1)

cmod.id = tmp_mid

print('\n{0:d}\n'.format(len(cmod.flux_bounds)) )
cmod.changeAllFluxBoundsWithValue(inf_bound, 'Infinity')
cmod.changeAllFluxBoundsWithValue(-inf_bound, '-Infinity')
cmod.changeAllFluxBoundsWithValue(numpy.inf, 'Infinity')
cmod.changeAllFluxBoundsWithValue(-numpy.inf, '-Infinity')

print('\n{0:d}\n'.format(len(cmod.flux_bounds)) )
slv.analyzeModel(cmod, lpFname=os.path.join(work_dir2, cmod.id), oldlpgen=False)

cmod.id = cmod.id.replace('_cplex','')+'.noinf'
print('\n{0:d}\n'.format(len(cmod.flux_bounds)) )
cmod.deleteAllFluxBoundsWithValue('Infinity')
cmod.deleteAllFluxBoundsWithValue('-Infinity')
print('\n{0:d}\n'.format(len(cmod.flux_bounds)) )

# We trust in SymPy to get the simplest representation of the rational by passing it a string!!!
# this code bypasses CBMPy's input checking and casting ...
for r_ in cmod.reactions:
    for rr_ in r_.reagents:
        rr_.coefficient = str(rr_.coefficient)
for fb_ in cmod.flux_bounds:
    fb_.value = str(fb_.value)

cmod.buildStoichMatrix(matrix_type='sympy')
CBWrite.exportModel(cmod, fmt='hformat', work_dir=work_dir2, use_rational=True)
