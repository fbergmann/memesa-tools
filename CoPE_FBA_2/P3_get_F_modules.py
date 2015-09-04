"""
Get F-modules
=============

This script determines the F-modules for a given model. We use a Python implementation of FluxModules (written by Arne Muller, now Arne Reimers). FluxModules requires as input a list of variable reactions. We use FVA output as input. Ideally, rational FVA is done, but this generally takes too long. 

This script accepts command line arguments (-m *model_name*  -l *sbml_level* -s *solver*  -c *ncores* -t *tolerance*) but you can also specify these arguments in this script

User input:
- *model_name* (string)
- *sbml_level* (integer)

Optional input:
- *solver* (string) [default = 'CPLEX'] Alternative is 'GLPK'
- *ncores* (integer) [default = 4]
- *tolerance* (float) [default = 10**-6, other values can result in incorrect determinition of the modules]

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

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 27, 2015
"""

from __future__ import division, print_function, absolute_import
import os, time, numpy,sys,getopt
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

### Set userdata here ####

model_name = 'toy_model'
sbml_level = 3
solver = 'CPLEX'
ncores = 4
tolerance = 10**-6

try:
    myopts, args = getopt.getopt(sys.argv[1:],"m:l:t:c:s:")
except getopt.GetoptError as e:
    print(str(e))
    print("Usage: %s -m model_name -l sbml_level -t tolerance -c ncores -s solver" % sys.argv[0])
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
    elif opt == '-t':
        tolerance = float(arg)
    elif opt == '-c':
        ncores = int(arg)
    elif opt == '-s':
        solver = arg
    else:
        print("Usage: %s -m model_name -l sbml_level -t tolerance -c ncores -s solver" % sys.argv[0])
        
if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')

from pyscescbm.CBVersion import __DEBUG__, __version__
import pyscescbm.fluxmodules.fluxmodules as fmod
from pyscescbm import CBRead, CBWrite, CBTools,CBMultiCore
from pyscescbm import CBSolver as slv

model_file = model_name + '.xml'
data_dir = os.path.join(cDir, 'data', model_name)
model_dir = os.path.join(data_dir,'models','sbml')
work_dir2 = os.path.join(data_dir,'models','h-format')
modules_dir = os.path.join(data_dir,'flux_modules')

def writeModule2CSV(flux_modules):
    """
    write calculate module to a csv file
    """
    file_path = os.path.join(modules_dir,"modules.txt")
    file_out = open(file_path,"w")        
    for module in flux_modules:
        if len(module) > 1:       # modules should contain more than one reaction
            for item in module[:-1]:
               file_out.write("{0}, ".format(item))
            file_out.write("{0}\n".format(module[-1]))
    file_out.close()

if not os.path.exists(work_dir2):
    os.makedirs(work_dir2)

if sbml_level == 3:
    print('Trying level 3')
    cmod = CBRead.readSBML3FBC(model_file, work_dir=model_dir)
else:
    cmod = CBRead.readSBML2FBA(model_file, work_dir=model_dir)
    
if solver.upper() == "CPLEX":
    fva_dat, fva_names = CBMultiCore.runMultiCoreFVA(cmod,procs=ncores)       
elif solver.upper() == "GLPK":
    fva_dat,fva_names = slv.glpk_FluxVariabilityAnalysis(cmod)    
else: 
    fva_dat, fva_names = CBMultiCore.runMultiCoreFVA(cmod,procs=ncores)     
    
lst_flux_modules = fmod.computeModules(cmod,tol=tolerance)    
writeModule2CSV(lst_flux_modules)
