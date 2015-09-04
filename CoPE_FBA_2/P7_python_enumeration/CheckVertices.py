from __future__ import division, print_function, absolute_import

input_dir = 'input'
temp_dir = 'temp'
output_dir = 'output'

import os
import numpy as np


try:
    import sympy
except ImportError:
    print("Error: SymPy is not installed")
    sys.exit()

Arr_b = []
Arr_A = []
Arr_Ab = []

def Negate(x):
    if x:
        x = -x
    return(x)
    

javapath = ":./jlinalg_0.5.jar" # Windows?
os.environ['CLASSPATH'] = javapath
print('CLASSPATH is now', os.environ['CLASSPATH'])
    
os.system("java CheckVertices {0:s} {1:d} {2:d} {3:s} | tee {4:s}.report.ver".format(os.path.join(temp_dir,"temple2"),25,8,os.path.join(temp_dir,"toy_model.2_split.noinf_r.ine.unpol.tidy.restored.onlyver"),os.path.join(temp_dir,"toy_model.2_split.noinf_r.ine"))) # TODO: Python implementation    
    
    
print("***********************")    
#sys.exit()    

### Import temple2
### temple 2 contains the Ax >= b data 
with open(os.path.join(temp_dir,"temple2"),'r') as infile:
    for line in infile:
       sline = line.strip().split(" ")
       sline_sympy = map(sympy.Rational,sline)        
       Arr_b.append(Negate(sline_sympy[0]))
       Arr_A.append(sline_sympy[1:])
       Arr_Ab.append(sline_sympy)

Arr_b = np.array(Arr_b)
Arr_A = np.array(Arr_A)
Arr_Ab = np.array(Arr_Ab)

rows = Arr_A.shape[0]
cols = Arr_A.shape[1]

#### Import vertices
with open(os.path.join(temp_dir,"toy_model.2_split.noinf_r.ine.unpol.tidy.restored.onlyver"),'r') as infile:
#with open(os.path.join(temp_dir,"Ecoli_iAF1260_ox.glc.1_split.noinf_r.ine.unpol.tidy.restored.onlyver"),'r') as infile:    
    for vertex in infile:        # for each vertex in infile    
       svertex = vertex.strip().split(" ")
       svertex_sympy = map(sympy.Rational,svertex)    
       Arr_v = np.transpose(svertex_sympy) 
       
       sol = np.dot(Arr_A,Arr_v) # computing Av
       
       ntight=0
       nloose=0
       fake_vertex=False
       Arr_sub_Ab = []           # submatrix that contains the set of tight inequalities constraints
       i=0
       ### Compare the outcome of Av to b
       for x,y in zip(sol,Arr_b):          
          if x==y:               # tight
              ntight+=1              
              Arr_sub_Ab.append(Arr_Ab[i])
          elif x>y:              # loose
              nloose+=1              
          else:                  # not in the polytope
              print("* IMPORTANT: This 'vertex' was not even in the polytope!")
              print("* This is a FAKE VERTEX")               
              fake_vertex=True        
              break
          i+=1

       if not fake_vertex:  
           ### Determine the rank of the submatrix which must be full rank
           print("* So {0:d} of the inequalities were tight. Checking the rank of this set...".format(ntight))
           Arr_sub_Ab = np.array(Arr_sub_Ab)
           u_, s_, v_ = np.linalg.svd(Arr_sub_Ab)
           matrix_rank = np.sum(s_ > 1e-10)
           print("* The submatrix has rank {0:d}".format(matrix_rank))           
           if matrix_rank == cols:  # full rank
               print("* This is a REAL VERTEX")           
           else: 
               print("* This is a FAKE VERTEX")
