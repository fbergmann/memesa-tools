"""
Rays statistics
===============

This script shows dictionaries of all rays which are not formed by splitting a reversible into a separate forward_suffix and backward_suffix reaction.

User input:
- *model_name* (string)
- *sbml_level* (integer)

Optional input:
- *forward_suffix* (string) [default = '_fwd']
- *backward_suffix* (string) [default = '_rev']

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

Written by Timo R. Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@vu.nl
Last Change: February 25, 2015
"""

from __future__ import division, print_function, absolute_import

import os,re,sys,h5py,numpy as np,sys,getopt
import pyscescbm as cbm

### DEFINE USER INPUT ###

model_name = 'toy_model'
sbml_level = 3
forward_suffix = '_fwd'
backward_suffix = '_rev'

#########################

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

from scripts import tools
from scripts import extract_fixed_network as fixed_network

cDir = os.getcwd()
data_dir = os.path.join(cDir, 'data',  model_name)
H_format_dir = os.path.join(data_dir,'models_subnetwork','h-format')
sbml_dir = os.path.join(data_dir,'models_subnetwork','sbml')
vertex_dir = os.path.join(data_dir,'cope_fba','subnetworks','vertex')
analysis_dir = os.path.join(cDir,'analysis',model_name)

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)   
    os.makedirs(os.path.join(analysis_dir,'subnetworks'))     
    os.makedirs(os.path.join(analysis_dir,'vertices'))

k=1
for file_in in sorted(os.listdir(vertex_dir)):                
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
            file_path = os.path.join(vertex_dir,file_in)
            f = h5py.File(file_path,'r')
        except Exception as er:
            print(er)
            print("Info: Extract {0:s} from the 7-zip file".format(file_path))
            sys.exit()            
        
        d = f['data']
        try:
            r = d['rays']
            HAVE_RAYS = True
        except KeyError:
            HAVE_RAYS = False
        
        if HAVE_RAYS:            
            L_r_ids = getReactionIdentifiers(file_in.replace('.None.hdf5',''),H_format_dir,delimiter = ',')    
            for j in range(r.shape[0]):                
                L_nonzero_indices = np.nonzero(r[j])[0]                
                D_fluxes = {L_r_ids[i] : r[j][i] for i in L_nonzero_indices}               
                if len(D_fluxes) == 2:
                    keys = sorted(D_fluxes)
                    if keys[0].endswith(forward_suffix) and keys[1].endswith(backward_suffix):
                        if keys[0][:len(forward_suffix)] != keys[1][:len(backward_suffix)]:                        
                            print(k, D_fluxes)  
                            k+=1
                    elif keys[1].endswith(forward_suffix) and keys[0].endswith(backward_suffix):
                        if not keys[1][:-len(forward_suffix)] != keys[0][:-len(backward_suffix)]:
                            print(k, D_fluxes)
                            k+=1
                    else:
                        print(k, D_fluxes)
                        k+=1
                    
                else:
                    print(k, D_fluxes)
                    k+=1
                    
                    
                    
                
                
                
