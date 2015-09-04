"""
Find shortest and cheapest path through the network
===================================================

User input:
- *model_name* (string)
- *sbml_level* (integer)

Optional input:
- *USE_PROTEIN_COSTS* (bool) [default = False]
- *rountto* (integer) [default = 7] 
- *L_protein_costs* (list) default = ['max.costs.csv','min.costs.csv','avg.costs.csv']

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, F.J. Bruggeman, and B. Teusink

Written by BG Olivier, Amsterdam, The Netherlands
E-mail: b.g.olivier@vu.nl
Last Change: February 25, 2015
"""

from __future__ import division, print_function, absolute_import

####### Settings #######
model_name = 'toy_model'
sbml_level = 3
USE_PROTEIN_COSTS = False
roundto = 7
L_protein_costs = ['max.costs.csv','min.costs.csv','avg.costs.csv']
########################

import os, sys, time, numpy,copy,getopt

try:
    myopts, args = getopt.getopt(sys.argv[1:],"m:l:")
except getopt.GetoptError as e:
    print(str(e))
    print("Usage: %s -m model_name -l sbml_level" % sys.argv[0])
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
    else:
        print("Usage: %s -m model_name -l sbml_level" % sys.argv[0]) 

if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
from scripts import tools
import pyscescbm as cbm

data_dir = os.path.join(cDir, 'data', model_name)
model_dir = os.path.join(data_dir,'models','sbml') 
analysis_dir = os.path.join(cDir,'analysis',model_name)
model_name = "{0:s}.xml".format(model_name)

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)   
    os.makedirs(os.path.join(analysis_dir,'subnetworks'))     
    os.makedirs(os.path.join(analysis_dir,'vertices'))

D_output = {}

if sbml_level == 3:
    cmod = cbm.CBRead.readSBML3FBC(model_name, work_dir=model_dir)    
else:
    cmod = cbm.CBRead.readSBML2FBA(model_name, work_dir=model_dir)    

if cmod == None:
   raise NameError("Specify the correct model name")
   
cmod.buildStoichMatrix()

L_r_ids = cmod.getReactionIds()
lp = cbm.CBSolver.cplx_MinimizeSumOfAbsFluxes(cmod,selected_reactions = L_r_ids) 
    
L_minSumAbsFluxes = []
for r in cmod.reactions:
    L_minSumAbsFluxes.append((r.getPid(), round(r.getValue(),roundto)))
    
D_output['SumAbsFluxes'] = {}
D_output['SumAbsFluxes']['min_SumAbsFluxes_flux'] = L_minSumAbsFluxes
D_output['SumAbsFluxes']['min_SumAbsFluxes'] = round(lp,roundto)

D_protein_costs = None
if USE_PROTEIN_COSTS:    
    cost_dir = os.path.join(data_dir,'protein_costs')    
    for cost_str in L_protein_costs:        
        D_protein_costs = tools.getProteinCosts('{0:s}.{1:s}'.format(model_name.split('.')[0],cost_str),cost_dir) 
        lp = cbm.CBSolver.cplx_MinimizeSumOfAbsFluxes(cmod,objective_coefficients=D_protein_costs,selected_reactions = L_r_ids)        
        L_minCost = []
        for r in cmod.reactions:
            L_minCost.append((r.getPid(), round(r.getValue(),roundto)))     
        cost_name = cost_str.rstrip('.csv')
        D_output[cost_name] = {}
        D_output[cost_name]['min_cost_flux'] = L_minCost
        D_output[cost_name]['min_cost'] = round(lp,roundto)  

L_poolDat = cbm.CBSolver.cplx_MinimizeNumActiveFluxes(cmod, selected_reactions=L_r_ids, populate=(0.0, 10, 300))
for i in range(1,len(L_poolDat[0])):           
    L_poolDat[0][i] = [round(elem,roundto) for elem in L_poolDat[0][i]]
    
cbm.CBWrite.printFBASolution(cmod)
D_output['opt'] =  cmod.getActiveObjective().getValue()
D_output['min_active_flux_pool'] = L_poolDat[0]
D_output['min_length'] = L_poolDat[1]
        
file_out = open(os.path.join(analysis_dir, 'MILP_results.py'), 'w')
file_out.write('results = \\\n')
file_out.write(cbm.CBTools.pprint.pformat(D_output))
file_out.flush()
file_out.close()
