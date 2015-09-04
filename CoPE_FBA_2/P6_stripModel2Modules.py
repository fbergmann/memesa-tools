"""
Strip the model into modules
============================

We use rational FBA and FluxModules (developed by Arne Reimers) output to "strip" the SBML model into subnetworks essential for optimizing the objective function. This allows for fast enumeration of the complete set of vertices. Stripped models are created in the following formats: INE, SBML, and MATLAB. Be aware of numerical issues in the SBML files (due to rounding).

This script accepts command line arguments (-m *model_name*  -b *inf_bound* -l *sbml_level*  -e *export_non_split* -s *model_splitting* -a *cytosol_abbreviation*  but you can also specify these arguments in this script

User input:
- *model_name* (string)
- *inf_bound* (integer/float)
- *sbml_level* (integer)

Optional input:
- *export_non_split* (bool) [default=False] generates subnetworks models of the original model in INE and SBML format
- *model_splitting* (bool) [default = False] generates a split SBML model

We use FluxModules on the complete flux space including futile cycles. Hence, we can obtain different subnetworks than we can find with CoPE-FBA because CoPE-FBA detects these futile cycles as rays and are therefore not included in the enumerated subnetworks. In practice, FluxModules returns modules which can be (a combination of) rays, linealities, and CoPE-FBA subnetworks.  Typically, these modules are CoPE-FBA subnetworks with coupled rays/linealities but they can also exist of one or multiple rays/linealities alone. The latter ones are not detected as subnetworks via CoPE-FBA (but as e.g. rays).

We added a "module input" and "module output" reaction to each subnetwork model which fixes the required input-output relationship (if any) of the subnetwork. Note that we also added "dummy species" to guarantee that both reactions must be used to reach steady state.

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

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 25, 2015
"""

from __future__ import division, print_function, absolute_import
import sys,getopt,os

from scripts import tools
############## USER INPUT #################

model_name = 'toy_model'
inf_bound = 1000
sbml_level = 3
export_non_split = False
model_splitting = False

try:
    myopts, args = getopt.getopt(sys.argv[1:],"m:b:l:e:s:")
except getopt.GetoptError as e:
    print(str(e))
    print("Usage: %s -m model_name -b inf_bound -l sbml_level -e export_non_split -s model_splitting" % sys.argv[0])
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
    elif opt == '-t':
        tolerance = float(arg)
    elif opt == '-e':
        export_non_split = bool(arg)
    elif opt == '-s':
        model_splitting = bool(arg)
    else:
        print("Usage: %s -m model_name -b inf_bound -l sbml_level -e export_non_split -s model_splitting" % sys.argv[0])

##########################################

def addFBAobjective(model,reaction_id,objective='maximize'):
    """     
    add an objective to the LP problem    
    
    Input:
      - *model*: (pyscescbm object)
      - *reaction_id*: (string)
      - *objective*: (string) [default= 'maximize']
    """
    r = model.getReaction(reaction_id)
    objective_flux = (1,reaction_id)
    objF = cbm.CBModel.Objective('max'+reaction_id,objective)
    objF.createFluxObjectives([objective_flux])
    model.addObjective(objF, active=True)
    return model

def addSpecies(model, species):
    """
    Input:
     - *model* (python object)
     - *species* (dict)
    """
    skeys = sorted(species)
    for S in skeys:        
        if 'compartment' in species[S]:
            comp = species[S]['compartment']
        else:
            print("Error! No compartment is given for species:\t",species[S]['id'])
 
        try:
            sObj = cbm.CBModel.Species(species[S]['id'], boundary=species[S]['boundary'],charge = species[S]['charge'],
                                       name=species[S]['name'], compartment=comp )
        except Exception as er:
            print(er)
            sys.exit()
            
        model.addSpecies(sObj)
    
def ParseRationalFBA(model_name,model_dir):
    """
    Parse rational FBA output generated with QSO PT EX version 2.5.0
    
    Input:
     - *model_name*
     - *model_dir*
    """
    ### Start by importing column identifiers
    column_identifiers = '{0:s}.noinf_r.columns.txt'.format(model_name) 
    D_column_ids = {}   
    try:
        file_in = open(os.path.join(model_dir,column_identifiers),'r')
    except IOError as er:
        print(er)
        sys.exit()
        
    n=1
    for line in file_in:
        try:
            D_column_ids['x{0:d}'.format(n)] = line.strip().split(',')[1]            
            n+=1
        except Exception as er:
            print(er)
            print(line)
    
    ### Import rational FBA output    
    try:
        file_in = open(os.path.join(model_dir,'{0:s}.noinf_r.ine.sol'.format(model_name) ),'r')
    except IOError as er:
        print(er)
        sys.exit()

    D_fluxes = {}
    IsVARS = False # variables (here, predicted flux values)
    for line in file_in:
       if line.startswith('slack'):
           IsVARS = False
           break
       elif IsVARS:
           sline = line.strip().split(' = ')       
           r_id = D_column_ids[sline[0]]
           J_r = sympy.Rational(sline[1])
           D_fluxes[r_id] = J_r
       elif line.startswith('VARS:'):
           IsVARS = True    
    return D_fluxes   
    
def ParseFluxModules(file_name,file_dir):
    """
    Parse output generated by FluxModules (developed by Arne Reimers)
    
    Input:
     - *file_name* (string)
     - *file_dir* (string)
    """
    try:     
        file_in = open(os.path.join(file_dir,file_name),'r')
    except Exception as er:
        print(er)
        sys.exit()   

    FluxModules = [line.strip().split(', ') for line in file_in]   
    return FluxModules   
    
def writeEFMToolInput(fba, fname=None, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    Write an FBA-LP in polynomial H-Format file. This is an improved version of `WriteModelHFormatFBA()`
    which it replaces but is kept for backwards compatability.

     - *fba* a PySCeS-CBM FBA object
     - *fname* [default=None] the output filename, fba.getPid() if not defined
     - *work_dir* [default=None] the output directory
     - *use_rational* [default=false] use rational numbers in output (requires sympy)
     - *fullLP* [default=True] include the default objective function as a maximization target
     - *format* [default='%s'] the number format string
     - *infinity_replace* [default=None] if defined this is the abs(value) of +-<infinity>
    """
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')
    M = fba    
    LHS = M.N.array.copy()
    RHS = [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
   
    RHS += [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
        
    objIdx = M.activeObjIdx
    LP = LHS
    del LHS, RHS
    
    if work_dir != None:
        assert os.path.exists(work_dir), '\nJanee ...'
        fname = os.path.join(work_dir, fname)

    if fname == None:
        fname = M.getPid().replace('.xml', '')
    if not use_rational:
        fname += '.m'
    else:
        fname += '_r.m'
        
    F = file(fname, 'w')
    F.write('% EFMTool input: {0}\n'.format(fname.split('/')[-1]) )    
    F.write('stru.stoich = [\n')           # stoichiometry matrix
    strW = format+' '
    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if not use_rational:
                if LP[r,c] == 0.0 or LP[r,c] == -0.0:
                    LP[r,c] = 0.0
                F.write(strW % LP[r,c])
            else:                
                F.write('{0},'.format(sympy.Rational(format % LP[r,c])))
        F.write('\n')
    F.write('];\n')

    F.write('\nstru.reversibilities = [ ') # reversibilities
    for r in M.reactions:
        if r.reversible:
            F.write('1 ')
        else:
            F.write('0 ')
    F.write('];\n')    
    
    F.write('\nstru.reactionNames = { ')  # reaction names
    nreactions = len(M.reactions)
    for i,r in enumerate(M.reactions):
        if i == (nreactions -1):
            F.write("'{0:s}' ".format(r.id) )
        else:
            F.write("'{0:s}', ".format(r.id) )
    F.write('};\n')
    F.close()
    return fname  
    
### Main ###

cytosol_abbrev = 'cytosol'
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
    
import os,copy,sys, numpy as np
import pyscescbm as cbm

_HAVE_SYMPY_ = None
try:
    import sympy
    _HAVE_SYMPY_ = True
except ImportError:
    print('Rational IO not available')
    _HAVE_SYMPY_ = False
    
__DEBUG__ = cbm.CBConfig.__CBCONFIG__['DEBUG']

if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')      

model_file = '{0:s}.xml'.format(model_name)
data_dir = os.path.join(cDir, 'data', model_name)
model_dir = os.path.join(data_dir,'models') 
fluxmodules_dir = os.path.join(data_dir,'flux_modules') # http://sourceforge.net/projects/fluxmodules/
subnetwork_dir = os.path.join(data_dir,'models_subnetwork')

if not os.path.exists(subnetwork_dir):
    os.makedirs(subnetwork_dir)     
    os.makedirs(os.path.join(subnetwork_dir,'sbml'))
    os.makedirs(os.path.join(subnetwork_dir,'h-format'))
    os.makedirs(os.path.join(subnetwork_dir,'MATLAB'))

if sbml_level  == 3 :
    try:
        cmod = cbm.CBRead.readSBML3FBC(model_file, os.path.join(model_dir,'sbml'))
    except:        
        cmod = cbm.CBRead.readSBML2FBA(model_file, os.path.join(model_dir,'sbml'))   
else:  # Try SBML L2
    cmod = cbm.CBRead.readSBML2FBA(model_file, os.path.join(model_dir,'sbml'))
    
cmod.setModifiedDate()
cmod.addModelCreator('Timo','Maarleveld', 'CWI Amsterdam and VU University Amsterdam')
cmod_original = copy.deepcopy(cmod)    

if model_splitting:  
    cmod2 = cbm.CBTools.splitReversibleReactions(cmod)
    cbm.CBWrite.writeSBML2FBA(cmod2, "{0:s}_split.xml".format(model_name), directory=os.path.join(model_dir,'sbml'))

D_fluxes = ParseRationalFBA(model_name,os.path.join(model_dir,'h-format'))
L_FluxModules = ParseFluxModules('modules.txt',fluxmodules_dir)

L_variable_ids = [item for sublist in L_FluxModules for item in sublist]
for r_id in L_variable_ids:
    with tools.suppress_stdout(): cmod.deleteReactionAndBounds(r_id)

### Determine number of "real" reversible reactions in the fixed part ### 
n=0
m=0
for r in cmod.reactions:
    if r.reversible:
       n+=1
       if r.getLowerBound() == -inf_bound and r.getUpperBound() == inf_bound:
           m+=1 
           
print("# Fixed reactions: {0:d}\n# Reversible reactions: {1:d}\n# Reversible reactions with -inf,inf bounds: {2:d}".format(len(cmod.reactions),n,m))

cmod = copy.deepcopy(cmod_original)
nmodules = len(L_FluxModules)
for i in range(nmodules):
    L_model_r_ids = cmod.getReactionIds()        
    flux_module = L_FluxModules[i]
    D_module_stoichiometry = tools.getExactInputOutputRelationship(cmod,flux_module,D_fluxes,quiet=True)  # rational input/output relationship (interface)
    #print(D_module_stoichiometry)    
    if D_module_stoichiometry != {}: # if NaJa == d        
        addSpecies(cmod,{ 'M_dummy_c' : {'id' : 'M_dummy_c', 'charge' : 0,'compartment' : cytosol_abbrev ,'boundary' : False, 'name': 'dummy'}}) # May 12, 2014
        D_reactions = {'R_module_output':{'lower': 1, 'reversible': False, 'reagents': [(-1,'M_dummy_c')], 'id': 'R_module_output', 'upper': 1},
                       'R_module_input':{'lower': 1, 'reversible': False, 'reagents': [(1,'M_dummy_c')], 'id': 'R_module_input', 'upper': 1}}         
        cbm.CBModelTools.addReactions(cmod,D_reactions)                             
        L_objectives = cmod.getObjectiveIds()        
        for obj in L_objectives:
            cmod.deleteObjective(obj)            
        addFBAobjective(cmod,'R_module_output')
    
    for r_id in L_model_r_ids:
        if (r_id not in flux_module):
            with tools.suppress_stdout(): cmod.deleteReactionAndBounds(r_id) # delete all reactions that do not occur in a module
    
    cmod.deleteNonReactingSpecies(simulate = False)
    cmod.setId("{0:s}.{1:d}".format(model_name,i+1))    

    cmod.deleteBoundsForReactionId(cmod.getActiveObjective().fluxObjectives[0].reaction)
    cmod2 = cbm.CBTools.splitReversibleReactions(cmod)
    cmod2.setId(cmod.getId()+'_split')

    if export_non_split:
        cmod.deleteAllFluxBoundsWithValue(inf_bound)
        cmod.deleteAllFluxBoundsWithValue(-inf_bound)       
    cmod2.deleteAllFluxBoundsWithValue(inf_bound)
    cmod2.deleteAllFluxBoundsWithValue(-inf_bound)
        
    # Write and Import Model to make sure that the order stays correct (It's weird I Know)       
    if export_non_split:
        cbm.CBWrite.writeSBML2FBA(cmod, "{0:s}.{1:d}.xml".format(model_name,i+1), directory=os.path.join(subnetwork_dir,'sbml'))   
        cmod = cbm.CBRead.readSBML2FBA("{0:s}.{1:d}.xml".format(model_name,i+1), os.path.join(subnetwork_dir,'sbml'))    
        cmod.deleteAllFluxBoundsWithValue(np.inf)
        cmod.deleteAllFluxBoundsWithValue(-np.inf)  
          
    cbm.CBWrite.writeSBML2FBA(cmod2, "{0:s}.{1:d}_split.xml".format(model_name,i+1), directory=os.path.join(subnetwork_dir,'sbml'))   
    cmod2 = cbm.CBRead.readSBML2FBA("{0:s}.{1:d}_split.xml".format(model_name,i+1), os.path.join(subnetwork_dir,'sbml'))
    cmod2.deleteAllFluxBoundsWithValue(np.inf)
    cmod2.deleteAllFluxBoundsWithValue(-np.inf)
    
    if D_module_stoichiometry != {}: # if NaJa == d != 0
        ### FIX module stoichiometry              
        for s_id in D_module_stoichiometry:
            s_stoich = D_module_stoichiometry[s_id]                
            if s_stoich < 0: 
                if export_non_split:
                    cmod.createReactionReagent('R_module_input',s_id,-s_stoich)               
                cmod2.createReactionReagent('R_module_input',s_id,-s_stoich)
            else:                
                if export_non_split:
                    cmod.createReactionReagent('R_module_output',s_id,-s_stoich)
                cmod2.createReactionReagent('R_module_output',s_id,-s_stoich)
            
    # We trust in SymPy to get the simplest representation of the rational by passing it a string!!!
    # this code bypasses CBMPy's input checking and casting ...
    cbm.CBWrite.writeSBML2FBA(cmod2, "{0:s}.{1:d}_split.correct.xml".format(model_name,i+1), directory=os.path.join(subnetwork_dir,'sbml'))    
    
    for r_ in cmod2.reactions:
        for rr_ in r_.reagents:
            rr_.coefficient = str(rr_.coefficient)        
    for fb_ in cmod2.flux_bounds:
        fb_.value = str(fb_.value)    

    cmod2.buildStoichMatrix(matrix_type='sympy')             
    cbm.CBWrite.writeModelHFormatFBA2(cmod2,fname = "{0:s}.{1:d}_split.noinf".format(model_name,i+1), work_dir = os.path.join(subnetwork_dir,'h-format'), use_rational=True)      
        
    if export_non_split:            
        cbm.CBWrite.writeSBML2FBA(cmod, "{0:s}.{1:d}.correct.xml".format(model_name,i+1), directory=os.path.join(subnetwork_dir,'sbml'))    
    
        for r_ in cmod.reactions:
            for rr_ in r_.reagents:
                rr_.coefficient = str(rr_.coefficient)  

        for fb_ in cmod.flux_bounds:
            fb_.value = str(fb_.value)   
  
        cmod.buildStoichMatrix(matrix_type='sympy')     
        cbm.CBWrite.writeModelHFormatFBA2(cmod, fname = "{0:s}.{1:d}.noinf".format(model_name,i+1), work_dir = os.path.join(subnetwork_dir,'h-format'), use_rational=True)        
    writeEFMToolInput(cmod2, fname = ("{0:s}.{1:d}_split.noinf".format(model_name,i+1)).replace('.','_'), work_dir = os.path.join(subnetwork_dir,'MATLAB'), use_rational=True)        
    cmod = copy.deepcopy(cmod_original) # reset cmod object
