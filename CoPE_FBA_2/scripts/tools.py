from __future__ import division, print_function, absolute_import
from contextlib import contextmanager
import os,sys

HAVE_SYMPY = False
try:
    import sympy
    HAVE_SYMPY = True
except ImportError:
    HAVE_SYMPY = False   

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def getProteinCosts(file_name,file_dir,IsHeader=True,delimiter=','):
    """
    Obtain the protein costs
    
    Input:
     - *file_name* (string)
    """
    D_protein_costs = None
    proteinCostsFile = os.path.join(file_dir,file_name)
    try:
        file_in = open(proteinCostsFile,'r')  # parse protein cost file        
    except Exception as er:        
        print(er)
        print("*** Add the protein costs document or set the USE_PROTEIN_COSTS argument to False ***")
        sys.exit()

    D_protein_costs = {}
    for line in file_in:
        if not IsHeader:
            line = line.strip()
            split_line = line.split(delimiter)    
            r_id = split_line[0]
            proteinCost = float(split_line[-1])  # takes the last column (HARDCODED)
            D_protein_costs[r_id] = proteinCost  # dictionary with protein costs                
        elif IsHeader:
            IsHeader=False    
    return D_protein_costs 

def getExactInputOutputRelationship(cmod,selected_reactions,D_fluxes,isboundary=False,quiet=False):     
    """
    Determine rational input-output relationship for a given set of reactions and rational fluxes
    
    Input:
    - *cmod*
    - *selected_reactions* (list)
    - *D_fluxes* (dict)        
    - *isboundary* (bool) [default = False]
    """
    D_stoichiometry = {}
    for r_id in selected_reactions:
        try:
            J_r = D_fluxes[r_id]            
        except KeyError:
            if not quiet:
                print("Warning: Reaction {0:s} not found in the provided flux dictionary. Flux is probably zero and therefore also set to zero".format(r_id) )
            J_r = 0
        if J_r:        
            r = cmod.getReaction(r_id)    
            L_stoichiometries = r.getStoichiometry()            
            for j in range(len(L_stoichiometries)):
                stoichiometry = L_stoichiometries[j][0]  # species stoichiometry
                s_id = L_stoichiometries[j][1]           # species identifier                
                s = cmod.getSpecies(s_id)                # get species object
                if isboundary:
                    if s_id in D_stoichiometry:                      
                        D_stoichiometry[s_id] += sympy.Rational(str(stoichiometry)).limit_denominator(10**20)*J_r
                    else:
                        D_stoichiometry[s_id] = sympy.Rational(str(stoichiometry)).limit_denominator(10**20)*J_r
                else:
                    if not s.is_boundary:                # ignore boundary species             
                        if s_id in D_stoichiometry:                      
                            D_stoichiometry[s_id] += sympy.Rational(str(stoichiometry)).limit_denominator(10**20)*J_r
                        else:
                            D_stoichiometry[s_id] = sympy.Rational(str(stoichiometry)).limit_denominator(10**20)*J_r                
                  
    for s_id in list(D_stoichiometry):        
        net_effect = D_stoichiometry[s_id]
        if not net_effect:
            D_stoichiometry.pop(s_id)
        else:
            D_stoichiometry[s_id] = net_effect         
    return D_stoichiometry
    
    
    
def getApproximateInputOutputRelationship(cmod,selected_reactions,D_fluxes,isboundary=False,roundto=7):     
    """
    Determine rational input-output relationship for a given set of reactions and rational fluxes
    
    Input:
    - *cmod*
    - *selected_reactions* (list)
    - *D_fluxes* (dict)        
    - *isboundary* (bool) [default = False]
    - *roundto* (integer) [default = 7]
    """
    D_stoichiometry = {}
    for r_id in selected_reactions:
        try:
            J_r = D_fluxes[r_id]            
        except KeyError:
            print("Warning: Reaction {0:s} not found in the provided flux dictionary. Flux is probably zero and therefore also set to zero".format(r_id) )
            J_r = 0
        if J_r:        
            r = cmod.getReaction(r_id)    
            L_stoichiometries = r.getStoichiometry()            
            for j in range(len(L_stoichiometries)):
                stoichiometry = L_stoichiometries[j][0]  # species stoichiometry
                s_id = L_stoichiometries[j][1]           # species identifier                
                s = cmod.getSpecies(s_id)                # get species object                
                if isboundary:
                    if s_id in D_stoichiometry:       
                        
                                       
                        D_stoichiometry[s_id] += stoichiometry*J_r
                    else:
                        D_stoichiometry[s_id] =  stoichiometry*J_r
                else:
                    if not s.is_boundary:                # ignore boundary species             
                        if s_id in D_stoichiometry:                      
                            D_stoichiometry[s_id] +=  stoichiometry*J_r
                        else:
                            D_stoichiometry[s_id] =  stoichiometry*J_r          
                  
    for s_id in list(D_stoichiometry):        
        net_effect = int(D_stoichiometry[s_id]*10**roundto)
        if not net_effect:
            D_stoichiometry.pop(s_id)
        else:
            D_stoichiometry[s_id] = net_effect/10**roundto
    return D_stoichiometry    
