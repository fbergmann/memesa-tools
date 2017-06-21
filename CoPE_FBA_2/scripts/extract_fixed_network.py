"""

This script determines attributes (length, cost, sum of abs. fluxes) of the fixed network

Needs: 
- Rational FBA (FVA) results of entire model
- (Genome-scale model, SBML L2 format, and Arne modules for sanity checks)

NEW: added dummy species

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 27, 2015
"""

from __future__ import division, print_function, absolute_import
import os,sys,time, sympy
# bgoli 2017
import cbmpy as cbm
from . import tools

cDir = os.getcwd()    

def ColumnCrossCheck(cmod,model_name,H_format_dir,delimiter=','):
    """
    Input:
     - *cmod* (CBMPy model object)
     - *model_filename* (string)
     - *delimiter* (string)
    
    # use the *.columns.txt file from the Steven output to cross check flux order
    """
    mismatch = []
    try:
        f_in = file(os.path.join(H_format_dir, '{0:s}.noinf_r.columns.txt'.format(model_name) ), 'r')
    except Exception as er:
        print(H_format_dir,model_name)
        print(er)
        sys.exit()
    indx_lst = []
    for line in f_in:
        l = [i.strip() for i in line.split(delimiter)]
        if len(l) == 2:
            indx_lst.append(l[1])
    f_in.close()
    for r in range(len(indx_lst)):
        match = bool(indx_lst[r] == cmod.N.col[r])
        print('{0} == {1}: {2}'.format(indx_lst[r], cmod.N.col[r], match) )
        if not match:
            mismatch.append((indx_lst[r], cmod.N.col[r]))
    if len(mismatch) > 0:
        print('\nColumn mismatch ({0:d})\n'.format(len(mismatch)) )
        os.sys.exit(1)
    else:
        print('\nColumn check successful, have a nice day :-)')
        time.sleep(2)        

def ParseFluxModules(module_name,module_dir):
    """
    
    Input:
     - *module_name* (string)
     - *module_dir* (string)
    
    http://sourceforge.net/projects/fluxmodules/
    """   
    try:     
       file_in = open(os.path.join(module_dir,module_name),'r')
    except Exception as er:
       print(er)
       sys.exit()

    L_module_r_ids = []
    for line in file_in :
        L_r_ids = line.strip().split(',')
        for j,name in enumerate(L_r_ids): # strip reaction identifiers                     
            L_r_ids[j] = name.strip()         
        L_module_r_ids += L_r_ids  
    return L_module_r_ids
    
def ParseRationalFBA(model_name,model_dir):
    """
    Input:
     - *model_name* (string)
     - *model_dir* (string)
    """
    column_identifiers = '{0:s}.noinf_r.columns.txt'.format(model_name)
    D_column_ids = {}
    try:
        file_in = open(os.path.join(model_dir,column_identifiers),'r')
    except IOError as er:
        print(er)        
        sys.exit()
    
    n=1
    for line in file_in:
        if line != '\n':
            try:
                D_column_ids['x{0:d}'.format(n)] = line.strip().split(',')[1]
                n+=1
            except Exception as er:
                print(er)
                print(line)
        
    try:
        file_in = open(os.path.join(model_dir,'{0:s}.noinf_r.ine.sol'.format(model_name) ),'r')
    except IOError as er:
        print(er)
        sys.exit()

    D_fluxes = {}
    IsVariables = False
    for line in file_in:
       if line.startswith('slack'):
           IsVariables = False
           break
       elif IsVariables:
           sline = line.strip().split(' = ')       
           r_id = D_column_ids[sline[0]]
           J_r = sympy.Rational(sline[1])
           D_fluxes[r_id] = J_r
       elif line.startswith('VARS:'):
           IsVariables = True
    return D_fluxes
    
def getFixedStats(model_name,data_dir,D_protein_costs,USE_COLUMN_CROSSCHECK=True,IS_FVA=False):
    """
    Input:
     - *model_name* (string)
     - *USE_COLUMN_CROSSCHECK* (boolean) [default = True]
    """    
    model_filename = '{0:s}.xml'.format(model_name)
    sbml_dir = os.path.join(data_dir,'models','sbml') 
    vertex_dir = os.path.join(data_dir,'cope_fba','whole_model','vertex')
    fluxmodules_dir = os.path.join(data_dir,'flux_modules') 
    cost_dir = os.path.join(data_dir,'protein_costs') 
    H_format_dir = os.path.join(data_dir,'models','h-format')   
    try:
        cmod = cbm.CBRead.readSBML3FBC(model_filename,sbml_dir)  
    except:
        cmod = cbm.CBRead.readSBML2FBA(model_filename,sbml_dir)          
    lp = cbm.CBSolver.analyzeModel(cmod,return_lp_obj = True)       

    if USE_COLUMN_CROSSCHECK:                  
        ColumnCrossCheck(cmod,model_name,H_format_dir) # Perform column cross check to make sure that identifiers match

    L_arne_module_r_ids = ParseFluxModules('modules.txt',fluxmodules_dir)        
    nvariable_fluxes = len(L_arne_module_r_ids)
    fixed_cost=0    
    fixed_sumAbsFluxes=0
    fixed_nfluxes=0       
    L_fixed_r_ids = []
    L_variable_r_ids = []
    D_fluxes = ParseRationalFBA(model_name,H_format_dir)      
    for r_id in D_fluxes:
        if r_id not in L_arne_module_r_ids: # == fixed
            J_r = D_fluxes[r_id]                          
            if J_r != 0:      
                L_fixed_r_ids.append(r_id)
                fixed_nfluxes+=1 
                fixed_sumAbsFluxes+=abs(J_r) 
                if D_protein_costs:
                    fixed_cost+=abs(J_r)*D_protein_costs[r_id]               
        else:
           L_variable_r_ids.append(r_id)
   
    if IS_FVA:
        try:
            columns_in = open(os.path.join(H_format_dir,'{0:s}.noinf_r.columns.txt'.format(model_name)))
            fva_in = open(os.path.join(vertex_dir,'{0:s}.noinf_r.ine.opt.fva'.format(model_name)))
        except Exception as er:
            print(er)
            sys.exit() 
             
        nJv=0
        nJf=0
        nJnf=0
        fixed_cost=0
        fixed_sumAbsFluxes=0
        for col,fva in zip(columns_in,fva_in):           
           J_r = sympy.Rational(fva.split(':')[1])  
           r_id = col.split(',')[1].strip()                  
           if 'FIXED' in fva:      
               assert r_id not in L_arne_module_r_ids, "Error: Unexpected behavior, fixed reactions cannot be part of a 'Fluxmodule'!"
               if J_r != 0:
                  nJf+=1          
                  fixed_sumAbsFluxes+=abs(J_r)
                  if D_protein_costs:
                      fixed_cost+=abs(J_r)*D_protein_costs[r_id]
               else:
                  nJnf += 1
                      
           elif 'VARIABLE' in fva:               
               nJv+=1
        assert nvariable_fluxes == nJv, "Error: Unexpected behavior, number of variable fluxes must be equal to the number of fluxes detected in all FluxModules together!"
        assert fixed_nfluxes == nJf, "Error: Unexpected behavior, number of variable fluxes must be equal to the number of fluxes detected in all FluxModules together!" 
    
    variable_IO = tools.getApproximateInputOutputRelationship(cmod,L_variable_r_ids,D_fluxes,isboundary=True)
    overal_IO = tools.getApproximateInputOutputRelationship(cmod, sorted(D_fluxes),D_fluxes,isboundary=True)     
    return fixed_nfluxes,nvariable_fluxes,float(fixed_cost),float(fixed_sumAbsFluxes),variable_IO,overal_IO  
