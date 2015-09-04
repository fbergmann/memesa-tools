"""
Determine vertex statistics
===========================

This script regenerates vertices from fixed and variable parts and determines length, sum of abs. fluxes and cost for each vertex.
Variable parts originate from enumerating vertices for each subnetwork (determiend with FluxModules).

User input:
- *model_name* (string)
- *sbml_level* (integer)

Optional input:
- *USE_PROTEIN_COSTS* (bool) [default = False]
- *TEST_EFM* (bool) [default = False]
- *DO_PLOTTING* (bool) [default = False]
- *DO_FULL_ENUMEREATION* [default = False] Set *DO_FULL_ENUMERATION* to True to determine the length, sum of abs. fluxes, and cost of each vertex. This is necessary for an (accurate) scatterplot. For models with many vertices etting *DO_FULL_ENUMERATION* to True results in much longer calculation times. The default setting of *DO_FULL_ENUMERATION* is False, since its not necessary to obtain hisograms of length, sum of abs. fluxes and cost of vertices.
- *forward_suffix* (string) [default = '_fwd'] 
- *backward_suffix* (string) [default = '_rev']

NOTE: We provide protein costs for the non-split model, thus use the correct suffixes for split reversible reactions for correct protein cost calculation. CBMPy uses ..._fwd and ..._rev. Example: R1_fwd and R1_rev

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

### DEFINE USER INPUT ###

model_name = 'toy_model'
sbml_level = 3

# Optional user input
TEST_EFM = False                # Test if enumerated (subnetwork) vertices are also EFMs 
USE_PROTEIN_COSTS = False
cost_suffix = 'max.costs.csv'   # suffix of the protein cost file
DO_PLOTTING = False             # histogram plotting
DO_FULL_ENUMERATION = False     # determine characteristics of each individual vertex
forward_suffix = "_fwd"
backward_suffix = "_rev"

# FIXED user input 

maxPoints= 10000 # max number of vertices of which we store both the length, sum of abs. fluxes, and cost
sfactor = 50.0   # Rounding floats: We convert floats to integers by multiplying it with a constant. This is used to quickly create an approximate histogram of vertex cost and length
input_function = 'R_module_input'
objective_function = 'R_module_output'

########################

import os,sys,itertools,numpy as np,h5py,operator,re,getopt

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

def sumproduct(*lists):
    return sum(reduce(operator.mul, data) for data in list(zip(*lists)))

def getReactionIdentifiers(model_file,model_dir,delimiter = ','):
    """
    Obtain the reaction identifiers corresponding to the model filename
    
    Input:
     - *model_file* (string)
     - *model_dir* (string)
     - *delimiter* (string) [default = ',']
    """
    file_path = os.path.join(model_dir,'{0:s}.noinf_r.columns.txt'.format(model_file) )
    try:        
        file_in = open(file_path,'r')
    except Exception as er:
        print(er)
        sys.exit()
        
    L_r_ids = []
    for line in file_in:        
        sline = line.strip().split(delimiter)        
        
        if len(sline) == 2:           # ignore empty lines
            r_id = sline[1]
            L_r_ids.append(r_id)    
    return L_r_ids
    
def GetVertexAttribute(D_current_values,to_add):
    """
    This function can be used to determine vertex attributes. Vertex attributates can be e.g. lengths, costs, or sum of absolute fluxes.
    
    Vertex attributes are determined via the determined subnetworks. In each subnetwork, different lengths/costs can be found. 
    
    We will use pathway length as an example to illustrate this function:
    - pathway lengths in subnetwork 1: [3,4]
    - pathway lengths in subnetwork 2: [2,2] 
    - pathway lengths in subnetwork 3: [2,3,4]
    - fixed network length = 45
 
    This results in 12 vertices with different lengths of which 52 (7+45) is the shortest and 55 (10+45) is the longest variable. Full enumeration is possible, but this becomes complicated if we have more subnetworks and/or subnetworks with more alternative pathways.
    
    Hence, we can use this function to iteratively add new subnetworks and count for each length the number of vertices. We initially start with a dictionary that represents the counts for the first subnetwork, *D_current_values* = {3:1,4:1}.
    We subsequently add the second subnetwork, *to_add* = [2,2]. This function returns the new lengths and counts, *D_updated_values* = {5:2,6:2}. 
    Then, we provide *D_updates_values* as input for the next addition where *to_add* = [2,3,4]. This gives *D_current_values* = {7:2,8:4,9:4,10:2}. 
    
    Input:
     - *D_current_values* (dict)
     - *to_add* (list)
    """
    minV = min(list(D_current_values))
    maxV = max(list(D_current_values))
    attr_values = list(range(minV,maxV+1))
    
    addMinV = min(to_add)
    addMaxV = max(to_add)

    D_updated_values = {}
    for l in range(min(attr_values)+addMinV,max(attr_values)+addMaxV+1):    
        total = 0
        i = min(attr_values)
        while i <= l - addMinV: 
            if i in D_current_values:
                total += to_add.count(l-i)*D_current_values[i]
            i+=1        
        if total:
            D_updated_values[l] = total
    return D_updated_values  
    
    
def writeAttributeData2CSV(module_number,r_ids,minL,maxL,minJ,maxJ,minC,maxC):
    """
    write attribute data to a CSV file
    """
    file_out = open(os.path.join(analysis_dir,'subnetworks',"attributes_{0:d}.csv".format(module_number)),"w")        
    file_out.write("reaction ids")
    file_out.writelines([",{0:s}".format(item)  for item in r_ids])          
    file_out.write('\n')
    for vertex in minL:
        file_out.write('minL')
        file_out.writelines([",{0}".format(item)  for item in vertex])   
        file_out.write('\n')
    for vertex in maxL:
        file_out.write('maxL')
        file_out.writelines([",{0}".format(item)  for item in vertex])   
        file_out.write('\n')
    for vertex in minJ:
        file_out.write('minJ')
        file_out.writelines([",{0}".format(item)  for item in vertex])   
        file_out.write('\n')
    for vertex in maxJ:
        file_out.write('maxJ')
        file_out.writelines([",{0}".format(item)  for item in vertex])   
        file_out.write('\n')
    if minC and maxC:
        for vertex in minC:
            file_out.write('minC')
            file_out.writelines([",{0}".format(item)  for item in vertex])   
            file_out.write('\n')
        for vertex in maxC:
            file_out.write('maxC')
            file_out.writelines([",{0}".format(item)  for item in vertex])   
            file_out.write('\n')
    file_out.close()      
        
def writeDistinctiveReactions2CSV(module_number,r_ids,attr,above_mean,below_mean,roundto = 7):
    """
    write reactions that cause the difference between both (or multiple) clusters to a file        
    
    Input:
     - *module_number* (int)
     - *r_ids* (list)
     - *attr* (string)
     - *above_mean* (list)
     - *below_mean* (list)
     - *roundto* (integer)
    """    
    absdiff = abs(np.mean(above_mean,0)-np.mean(below_mean,0))
    diff = np.mean(above_mean,0)-np.mean(below_mean,0)        
    reaction_influences = sorted(zip(diff,r_ids), key=operator.itemgetter(0),reverse=True)    
    diff_sorted = [list(x) for x in zip(*sorted(reaction_influences))]
    filename = os.path.join(analysis_dir,'subnetworks','{0:d}_{1:s}_explaining_reactions.csv'.format(module_number,attr))
    file_out = open(filename,'w')       
    for i in range(len(r_ids)):        
        if abs(diff_sorted[0][i]) > 10**-roundto: # ignore values smaller than the precision of the simulator
            file_out.write("{0:s},{1}\n".format(diff_sorted[1][i],diff_sorted[0][i]))   
    file_out.close()    
    
def IsEFM(N_matrix,v):
    """ Use the rank test to determine if a enumerated vertex is also an EFM """
    N_submatrix = np.delete(N_matrix,np.where(v == 0)[0],1)    
    ncolumns = N_submatrix.shape[1]  
    u_, s_, v_ = np.linalg.svd(N_submatrix)
    rank = np.sum(s_ > 1e-10)
    
    nullity = ncolumns - rank    
    IsEFM = False
    if nullity == 1:
        IsEFM = True
    else:
        print("Warning: No EFM. Rank = {0:d}. Nullity is {0:d}.".format(rank,nullity) )
    return IsEFM
    
def getVariableNetworkStats():
    print('-------------------------------------------------------------')
    L_variable_lengths = []
    L_variable_costs = []
    L_variable_sumAbsFluxes = []
    D_IO = {}
    minSumAbsFluxes = fixed_sumAbsFluxes
    maxSumAbsFluxes = fixed_sumAbsFluxes
    minCost = fixed_cost
    maxCost = fixed_cost
    minLength = fixed_length
    maxLength = fixed_length
    nvertices_network = 1
    nMinLength = 1
    nMinCost = 1
    nMinSumAbsFluxes = 1
    n=1
    if USE_PROTEIN_COSTS:
       file_out = open(os.path.join(analysis_dir,'subnetworks','subnetwork_stats_{0:s}'.format(cost_suffix) ),'w')
       file_out.write("Module,dL,dJ,dC\n")
    else:
       file_out = open(os.path.join(analysis_dir,'subnetworks','subnetwork_stats.csv'),'w')
       file_out.write("Module,dL,dJ\n")
    for file_in in sorted(os.listdir(work_dir)):                
        if file_in.endswith('hdf5') and model_name in file_in:                
            regex = re.findall('\.\d+[\_\.]',file_in) # \_            
            module_number = int(regex[0][1:-1])
            with tools.suppress_stdout():     
                if sbml_level == 3:
                    try:
                        cmod = cbm.CBRead.readSBML3FBC(file_in.replace('.None.hdf5','.correct.xml'),sbml_dir)     
                    except: 
                        cmod = cbm.CBRead.readSBML2FBA(file_in.replace('.None.hdf5','.correct.xml'),sbml_dir)                      
                else:
                    cmod = cbm.CBRead.readSBML2FBA(file_in.replace('.None.hdf5','.correct.xml'),sbml_dir)                  
                     
            D_attributes = {'length':[],'cost':[],'sumAbsFluxes':[]}
            try:
                file_path = os.path.join(work_dir,file_in)
                f = h5py.File(file_path,'r')
            except Exception as er:
                print(er)
                print("Info: Extract {0:s} from the 7-zip file".format(file_path))
                sys.exit()            
            
            d = f['data']
            try:
                v = d['vertices']
                HAVE_VERTICES = True
            except KeyError:
                HAVE_VERTICES = False
            
            if HAVE_VERTICES:
                cmod.buildStoichMatrix()    
                L_r_ids = getReactionIdentifiers(file_in.replace('.None.hdf5',''),H_format_dir,delimiter = ',')                  
                if TEST_EFM:
                    N_matrix = np.ones([len(cmod.N.row),len(cmod.N.col)])
                    for i,r_id in enumerate(cmod.N.col):
                        j = L_r_ids.index(r_id)            
                        N_matrix[:,j] = cmod.N.array[:,i]           
                         
                objective_index = L_r_ids.index(objective_function)    
                input_index = L_r_ids.index(input_function)        
                L_r_ids = np.delete(L_r_ids,[input_index,objective_index])                   
                
                if USE_PROTEIN_COSTS:
                    Arr_protein_costs = np.ones(len(L_r_ids))
                    for i,r_id in enumerate(L_r_ids):  # hardcoded names for forward and backward reaction ...         
                        r_id = r_id.replace(backward_suffix,'')
                        r_id = r_id.replace(forward_suffix,'')
                        Arr_protein_costs[i] = D_protein_costs[r_id]        
                
                nvertices_network*=v.shape[0]
                L_IO_relationships = []
                minC = 999999
                maxC = 0
                minL = 999999
                maxL = 0
                minJ = 999999
                maxJ = 0
                L_maxC = None
                L_minC = None
                for i in range(v.shape[0]):
                    assert v[i,objective_index]/v[i,input_index] == 1, "Error: yield must be 1!" ### DOUBLE CHECK IF YIELD = 1 ###
                    Arr_vertex = np.delete(v[i], [input_index,objective_index])
                    L_nonzero_indices = np.nonzero(Arr_vertex)[0]                
                    
                    if TEST_EFM:
                    
                        IsEFM(N_matrix,v[i]) 
                                
                    D_fluxes = {L_r_ids[i] : Arr_vertex[i] for i in L_nonzero_indices}
                    vertex_IO = tools.getApproximateInputOutputRelationship(cmod,sorted(D_fluxes),D_fluxes) 
                    L_IO_relationships.append(vertex_IO)     
                    
                    ### Attributes ###
                    vertex_len = len(L_nonzero_indices)
                    D_attributes['length'].append(vertex_len)                               
                   
                    vertex_sumAbsFlux = np.abs(Arr_vertex).sum() # we use abs if we analyze a non-split model with reversible reactions
                    D_attributes['sumAbsFluxes'].append(vertex_sumAbsFlux)
                    
                    ### min-max vertex Length ###
                    if vertex_len < minL:
                        L_minL = [Arr_vertex]
                        minL = vertex_len
                    elif vertex_len == minL:
                        L_minL.append(Arr_vertex)
                    if vertex_len > maxL:
                        L_maxL = [Arr_vertex]
                        maxL = vertex_len
                    elif vertex_len == maxL:
                        L_maxL.append(Arr_vertex)               
                    
                    ### min-max vertex sum of abs. fluxes ###
                    if vertex_sumAbsFlux < minJ:
                        L_minJ = [Arr_vertex]
                        minJ = vertex_sumAbsFlux
                    elif vertex_sumAbsFlux == minJ:
                        L_minJ.append(Arr_vertex)
                    if vertex_sumAbsFlux > maxJ:
                        L_maxJ = [Arr_vertex]
                        maxJ = vertex_sumAbsFlux
                    elif vertex_sumAbsFlux == maxJ:
                        L_maxJ.append(Arr_vertex)       
                    
                    ### min-max vertex cost ###
                    if USE_PROTEIN_COSTS:            
                        Arr_vertex_costs = np.multiply(np.abs(Arr_vertex),Arr_protein_costs)   # piecewise multiplication 24 july
                        vertex_cost = Arr_vertex_costs.sum()                       
                        if vertex_cost < minC:
                            L_minC = [Arr_vertex]
                            minC = vertex_cost
                        elif vertex_cost == minC:
                            L_minC.append(Arr_vertex)
                        if vertex_cost > maxC:
                            L_maxC = [Arr_vertex]
                            maxC = vertex_cost
                        elif vertex_cost == maxC:
                            L_maxC.append(Arr_vertex)
                        D_attributes['cost'].append(vertex_cost)            
                
                writeAttributeData2CSV(module_number,L_r_ids,L_minL,L_maxL,L_minJ,L_maxJ,L_minC,L_maxC)
                #print(L_IO_relationships)
                D_module_IO = {k: v for k,v in  L_IO_relationships[0].items()}

                #print(L_IO_relationships)
                for i in range(1,len(L_IO_relationships)):
                    assert sorted(D_module_IO) == sorted(L_IO_relationships[i]),"Error: ... : \n{0}\n{1}".format(D_module_IO,L_IO_relationships[i])
                    for s_id in D_module_IO:
                        assert str(D_module_IO[s_id]) == str(L_IO_relationships[i][s_id]), "Error! Input-output relationship does not match: \n{0}\n{1}".format(D_module_IO[s_id],L_IO_relationships[i][s_id])
                                        
                    #assert L_IO_relationships[i] == D_module_IO, "Error! Input-output relationship does not match: \n{0}\n{1}".format(} ,D_module_IO)
               
                for s_id in list(D_module_IO):
                    if s_id not in D_IO:
                        D_IO[s_id] = D_module_IO[s_id]
                    else:
                        D_IO[s_id] += D_module_IO[s_id]    
                
                print("#  Module {0:s} with {1:d} biomass vertices; total vertices {2:d}".format(file_in.replace('.None.hdf5',''), v.shape[0], nvertices_network ))
                module_maxL = max(D_attributes['length'])
                module_minL = min(D_attributes['length'])
                print("   Length: {0:d}-{1:d}".format(module_minL,module_maxL) )
                dL = module_maxL - module_minL
                minLength += module_minL 
                maxLength += module_maxL 
                nMinLength *= D_attributes['length'].count(module_minL )
                L_variable_lengths.append(D_attributes['length'])              
                #if dL:
                #    plt.figure()        
                #    plt.hist(D_attributes['length'],bins = dL)
                
                module_maxSAF = max(D_attributes['sumAbsFluxes'])
                module_minSAF = min(D_attributes['sumAbsFluxes'])
                print("   SumAbsFluxes: {0}-{1}".format (module_minSAF,module_maxSAF) )        
                dSAF = module_maxSAF-module_minSAF
                minSumAbsFluxes += module_minSAF
                maxSumAbsFluxes += module_maxSAF
                
                str_sumAbsFluxes = [str(J) for J in D_attributes['sumAbsFluxes']]
                nMinSumAbsFluxes *= str_sumAbsFluxes.count(str(module_minSAF))

                L_variable_sumAbsFluxes.append(D_attributes['sumAbsFluxes']) 
                #if dSAF > 1:        
                #    plt.figure()        
                #    plt.hist(D_attributes['sumAbsFluxes'],bins = dSAF)      
                
                meanSumAbsFluxes = np.mean(D_attributes['sumAbsFluxes'])        
                above_mean = []
                below_mean = []
                for i in range(v.shape[0]):      
                    Arr_vertex = np.delete(v[i], [input_index,objective_index])
                    vertexSumAbsFluxes = np.abs(Arr_vertex).sum()
                    if vertexSumAbsFluxes > meanSumAbsFluxes:
                        above_mean.append(Arr_vertex)
                    else:
                        below_mean.append(Arr_vertex)         
                        
                if module_maxSAF - module_minSAF > 10**-6:                
                    writeDistinctiveReactions2CSV(module_number,L_r_ids,"sumAbsFluxes",above_mean,below_mean) 
                    
                if USE_PROTEIN_COSTS:
                    module_maxC = max(D_attributes['cost'])
                    module_minC = min(D_attributes['cost'])
                    print("   Cost: {0}-{1}".format(module_minC,module_maxC) )
                    dC = module_maxC-module_minC
                    minCost += module_minC
                    maxCost += module_maxC
                    
                    str_cost = [str(J) for J in D_attributes['cost']]                
                    nMinCost *= str_cost.count(str(module_minC))                
                    L_variable_costs.append(D_attributes['cost'])            
                       
                    #if dC > 1:        
                    #  plt.figure()        
                    #  plt.hist(D_attributes['cost'],bins = dC)     
                             
                    meanCost = np.mean(D_attributes['cost'])                                   
                    above_mean = []
                    below_mean = []
                    for i in range(v.shape[0]):      
                        Arr_vertex = np.delete(v[i], [input_index,objective_index])
                        Arr_vertex_costs = np.multiply(np.abs(Arr_vertex),Arr_protein_costs)   # piecewise multiplication 
                        vertex_cost = Arr_vertex_costs.sum()                    
                        if vertex_cost > meanCost :
                            above_mean.append(Arr_vertex_costs)
                        else:
                            below_mean.append(Arr_vertex_costs)   
                    if module_maxC - module_minC > 10**-6:                      
                        writeDistinctiveReactions2CSV(module_number,L_r_ids,"cost",above_mean,below_mean)                                  
                    file_out.write("{0:s},{1:d},{2},{3}\n".format(file_in.replace('.None.hdf5',''),dL,dSAF,dC))              
                else:
                    file_out.write("{0:s},{1:d},{2}\n".format(file_in.replace('.None.hdf5',''),dL,dSAF))              
            n+=1
    file_out.close()
    print('\n-------------------------------------------------------------')
    print("Min-Max Length: {0:d} {1:d} {2}".format(minLength,maxLength,100*(maxLength-minLength)/float(minLength)))   
    print("Min-Max SumOfFluxes: {0} {1} {2}".format(minSumAbsFluxes,maxSumAbsFluxes,100*(maxSumAbsFluxes-minSumAbsFluxes)/minSumAbsFluxes))
    if USE_PROTEIN_COSTS:
        print("Min-Max Cost: {0} {1} {2}".format(minCost,maxCost,100*(maxCost-minCost)/minCost))
    print('\n-------------------------------------------------------------')     
    
    print("\n# vertices in min(length)={0:d} solution space: {1:d}".format(minLength,nMinLength))
    print("# vertices in min(sumAbsFluxes)={0} solution space: {1:d}".format(minSumAbsFluxes, nMinSumAbsFluxes))
    if USE_PROTEIN_COSTS:
        print("# vertices in min(cost)={0} solution space: {1:d}\n".format(minCost, nMinCost))
    #plt.show()
    D_stats = {}
    D_stats['minL'] = minLength
    D_stats['maxL'] = maxLength
    D_stats['minSAF'] = minSumAbsFluxes
    D_stats['maxSAF'] = maxSumAbsFluxes
    D_stats['minC'] = minCost
    D_stats['maxC'] = maxCost
    return D_IO,L_variable_lengths,L_variable_sumAbsFluxes, L_variable_costs,nvertices_network,D_stats


import pyscescbm as cbm
if DO_PLOTTING:
    import matplotlib as mpl, matplotlib.pyplot as plt
    import matplotlib.mlab as mlab

from scripts import tools
from scripts import extract_fixed_network as fixed_network

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

data_dir = os.path.join(cDir, 'data',  model_name)
H_format_dir = os.path.join(data_dir,'models_subnetwork','h-format')
sbml_dir = os.path.join(data_dir,'models_subnetwork','sbml')
work_dir = os.path.join(data_dir,'cope_fba','subnetworks','vertex')
analysis_dir = os.path.join(cDir,'analysis',model_name)

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)   
    os.makedirs(os.path.join(analysis_dir,'subnetworks'))     
    os.makedirs(os.path.join(analysis_dir,'vertices'))
    
if os.listdir(work_dir) == []:
    raise UserWarning("No files to analyze in directory {0:s}".format(work_dir))
    
### Get protein costs ###    
D_protein_costs = None
if USE_PROTEIN_COSTS:
    cost_dir = os.path.join(data_dir,'protein_costs')    
    D_protein_costs =  tools.getProteinCosts('{0:s}.{1:s}'.format(model_name.split('.')[0],cost_suffix),cost_dir)     


### Get fixed network stats ###
fixed_length = 0
fixed_sumAbsFluxes = 0
fixed_cost = 0
#with tools.suppress_stdout():         
fixed_length,variable_length,fixed_cost,fixed_sumAbsFluxes, variable_IO,overall_IO = fixed_network.getFixedStats(model_name,data_dir,D_protein_costs,USE_COLUMN_CROSSCHECK=True,IS_FVA=False)

print("# fixed reactions: {0:d}".format(fixed_length) )
D_input_output,L_variable_lengths,L_variable_sumAbsFluxes,L_variable_costs,nvertices_network,D_stats = getVariableNetworkStats()

### Check Input-output relationship###
to_remove = []
for s_id in D_input_output:    
    if D_input_output[s_id]: # some species can have an overal network stoichiometry of 0 (+1 in module X and -1 in module Y)        
        assert abs(D_input_output[s_id] - variable_IO[s_id]) < 10**-6, "Error! net stoichiometry mismatch for {0:s} : {1} {2}".format(s_id,D_input_output[s_id],variable_IO[s_id])     
    
print("\nNet stoichiometry\n==================")
for s_id in sorted(overall_IO):
    print("{0:s} {1}".format(s_id,overall_IO[s_id]))    
print('\n-------------------------------------------------------------\n')

### Lengths ###
D_current_values = {}
for length in L_variable_lengths[0]:
    if (fixed_length + length) not in D_current_values: # July 08, important modification!
        D_current_values[fixed_length + length] = 1
    else:
        D_current_values[fixed_length + length] += 1
for to_add in L_variable_lengths[1:]:    
    D_current_values = GetVertexAttribute(D_current_values,to_add)

L_lengths = sorted(D_current_values)
L_counts = []
file_path = os.path.join(analysis_dir,'vertices','hist_len_counts.csv')
file_out = open(file_path,'w' )
for length in L_lengths :
   count = D_current_values[length]   
   file_out.write("{0:d},{1:d}\n".format(length, count) )
   L_counts.append(count)
file_out.close()
print("Info: hist length counts data successfully saved at {0}".format(file_path) )

if DO_PLOTTING:
    plt.figure()
    n, bins, patches = plt.hist(L_lengths,weights = L_counts,bins = D_stats['maxL']-D_stats['minL'],align='left')
    plt.xlim((D_stats['minL'],D_stats['maxL']))
    mu = sumproduct(L_lengths,L_counts)/float(sum(L_counts))    
    var = sumproduct((np.array(L_lengths)-mu)**2,L_counts)/(sum(L_counts)-1.0)        
    y = mlab.normpdf(bins, mu, var**0.5)  # add normal distribution with sample mean and standard deviation
    plt.plot(bins, np.array(y)*sum(L_counts), 'r--') 
    
### Sum Absolute Fluxes ###
L_scaled_sumAbsFluxes1 = np.rint(np.array(L_variable_sumAbsFluxes[0])*sfactor).astype(int).tolist() # make integers of variable costs (by multiplying with a constant > 1 )
D_current_sumAbsFluxes = {}
for Js in L_scaled_sumAbsFluxes1:
    if Js not in D_current_sumAbsFluxes:
        D_current_sumAbsFluxes[Js] = 1 # we add fixed cost in a later phase
    else:
        D_current_sumAbsFluxes[Js] += 1 
for to_add in L_variable_sumAbsFluxes[1:]:
    L_scaled_sumAbsFluxes2 = np.rint(np.array(to_add)*sfactor).astype(int).tolist() # make integers of variable costs (by multiplying with a constant > 1 )
    D_current_sumAbsFluxes = GetVertexAttribute(D_current_sumAbsFluxes,L_scaled_sumAbsFluxes2)

file_path = os.path.join(analysis_dir,'vertices','hist_sumAbsFluxes_counts.csv')       
file_out = open(file_path,'w' )
L_var_sumAbsFluxes = sorted(D_current_sumAbsFluxes)
L_counts = []
for Js in L_var_sumAbsFluxes:
   L_counts.append(D_current_sumAbsFluxes[Js])
   file_out.write("{0},{1}\n".format(fixed_sumAbsFluxes + Js/sfactor, D_current_sumAbsFluxes[Js]) )
print("Info: hist sum of abs. fluxes data successfully saved at {0:s}".format(file_path) )
file_out.close()
 
if DO_PLOTTING:
    L_sumAbsFluxes = fixed_sumAbsFluxes + np.array(L_var_sumAbsFluxes)/sfactor
    plt.figure()
    plt.hist(L_sumAbsFluxes,weights = L_counts,bins = np.ceil(D_stats['maxSAF']-D_stats['minSAF']),align='left')
    plt.xlim((D_stats['minSAF'],D_stats['maxSAF']))
    
if USE_PROTEIN_COSTS:
    L_scaled_costs1 = np.rint(np.array(L_variable_costs[0])*sfactor).astype(int).tolist() # make integers of variable costs (by multiplying with a constant > 1 )
    D_current_costs = {}
    for cost in L_scaled_costs1:
        if cost not in D_current_costs:
            D_current_costs[cost] = 1 # we add fixed cost in a later phase
        else:
            D_current_costs[cost] += 1 
    for to_add in L_variable_costs[1:]:    
        L_scaled_costs2 = np.rint(np.array(to_add)*sfactor).astype(int).tolist() # make integers of variable costs (by multiplying with a constant > 1 )
        D_current_costs = GetVertexAttribute(D_current_costs,L_scaled_costs2)

    file_path = os.path.join(analysis_dir,'vertices','hist_cost_counts_{0:s}.csv'.format(cost_suffix.replace('.costs.csv','') ) )   
    file_out = open(file_path,'w' )
    L_var_costs = sorted(D_current_costs)
    L_counts = []
    for cost in L_var_costs:      
       L_counts.append(D_current_costs[cost])
       file_out.write("{0},{1}\n".format(fixed_cost + cost/sfactor, D_current_costs[cost]))
    file_out.close()
    print("Info: hist cost data successfully saved at {0:s}".format(file_path) )
    
    if DO_PLOTTING:
        L_costs = fixed_cost + np.array(L_var_costs)/sfactor
        plt.figure()
        plt.hist(L_costs,weights = L_counts,bins = np.ceil(D_stats['maxC']-D_stats['minC']),align = 'left')
        plt.xlim((D_stats['minC'],D_stats['maxC']))

if DO_FULL_ENUMERATION:
    Arr_vertex_lengths = np.zeros(nvertices_network,dtype=int)
    nMinLength = 0
    for i,pair in enumerate(itertools.product(*L_variable_lengths)):
        length = sum(pair) + fixed_length
        if length == D_stats['minL']:
            nMinLength += 1
        Arr_vertex_lengths[i] = int(length)
        if i and not i%10000000:
            print(i)
    print("# minLength: {0:d}".format(nMinLength) )
    from collections import Counter
    D_len_counts = Counter(Arr_vertex_lengths)    
    file_path = os.path.join(analysis_dir,'vertices','hist_len_counts_full.csv')
    file_out = open(file_path,'w' )
    for length in sorted(D_len_counts):
        file_out.write("{0:d},{1:d}\n".format(length, D_len_counts[length]) )
    file_out.close()
    print("Info: hist cost counts data successfully saved at {0:s}".format(file_path) )
    if DO_PLOTTING:
        plt.figure()
        n, bins, patches = plt.hist(Arr_vertex_lengths,bins = D_stats['maxL']-D_stats['minL'],align='left')        
        y = mlab.normpdf(bins, np.mean(Arr_vertex_lengths), np.std(Arr_vertex_lengths))  # add normal distribution with sample mean and standard deviation
        plt.plot(bins, np.array(y)*len(Arr_vertex_lengths), 'r--')    
    
    Arr_vertex_sumAbsFluxes = np.zeros(nvertices_network)
    nMinSumAbsFluxes = 0
    for i,pair in enumerate(itertools.product(*L_variable_sumAbsFluxes)):
        sumAbsFluxes = sum(pair) + fixed_sumAbsFluxes
        if str(sumAbsFluxes) == str(D_stats['minSAF']):
            nMinSumAbsFluxes += 1
        Arr_vertex_sumAbsFluxes[i] = sumAbsFluxes
        if i and not i%10000000:
            print(i)
    print("# minSumAbsfluxes: {0:d}".format(nMinSumAbsFluxes) )
    D_sumAbsFluxes_counts = Counter(np.floor(Arr_vertex_sumAbsFluxes))
    if USE_PROTEIN_COSTS:
        file_path = os.path.join(analysis_dir,'vertices','hist_sumAbsFluxes_counts_{0:s}_full.csv'.format(cost_suffix.replace('.costs.csv','') ) )
    else:
        file_path = os.path.join(analysis_dir,'vertices','hist_sumAbsFluxes_counts_full.csv')
    file_out = open(file_path,'w' )
    for SAF in sorted(D_sumAbsFluxes_counts):
        file_out.write("{0},{1}\n".format(SAF, D_sumAbsFluxes_counts[SAF]))
    file_out.close()
    print("Info: hist cost counts data successfully saved at {0:s}".format(file_path) )
    if DO_PLOTTING:
        plt.figure()       
        plt.hist(Arr_vertex_sumAbsFluxes,bins = np.ceil(D_stats['maxSAF']-D_stats['minSAF']),align='left')
        plt.xlim((D_stats['minSAF'],D_stats['maxSAF']))
        
    if USE_PROTEIN_COSTS:  
        Arr_vertex_costs = np.zeros(nvertices_network)
        nMinCost = 0
        for i,pair in enumerate(itertools.product(*L_variable_costs)):
            cost = sum(pair) + fixed_cost
            if str(cost) == str(D_stats['minC']):
                nMinCost += 1
            Arr_vertex_costs[i] = cost
            if i and not i%10000000:
                print(i)
        print("# minCost {0:d}:".format(nMinCost) )
        D_cost_counts = Counter(np.floor(Arr_vertex_costs))
        if USE_PROTEIN_COSTS:
            file_path = os.path.join(analysis_dir,'vertices','hist_cost_counts_{0:s}_full.csv'.format(cost_suffix.replace('.costs.csv','') ) )
        else:
            file_path = os.path.join(analysis_dir,'vertices','hist_cost_counts_full.csv')
        file_out = open(file_path,'w' )
        for cost in sorted(D_cost_counts):
            file_out.write("{0},{1}\n".format(cost, D_cost_counts[cost]) )
        file_out.close()
        print("Info: hist cost counts data successfully saved at {0:s}".format(file_path) )
        if DO_PLOTTING:
            plt.figure()
            plt.hist(Arr_vertex_costs,bins = np.ceil(D_stats['maxC']-D_stats['minC']),align='left')
            plt.xlim((D_stats['minC'],D_stats['maxC']))
        
if DO_PLOTTING:
    plt.show()

if DO_FULL_ENUMERATION:
    if USE_PROTEIN_COSTS:
        file_path = os.path.join(analysis_dir,'vertices','vertex_stats_{0:s}.csv'.format(cost_suffix.replace('.costs.csv','') ) )
    else:
        file_path = os.path.join(analysis_dir,'vertices','vertex_stats.csv')
    file_out = open(file_path,'w' )

    modulo = 1
    if nvertices_network > maxPoints:
        modulo = nvertices_network//maxPoints  # Use double // for Python 3
 
    for i in range(nvertices_network):        
        if not i%modulo:             
            if USE_PROTEIN_COSTS:
                file_out.write("{0},{1},{2}\n".format(Arr_vertex_lengths[i],Arr_vertex_sumAbsFluxes[i],Arr_vertex_costs[i]))
            else:
                file_out.write("{0},{1}\n".format(Arr_vertex_lengths[i],Arr_vertex_sumAbsFluxes[i]))                
    file_out.close()
    print("Info: vertex data successfully saved at {0:s}".format(file_path) )
