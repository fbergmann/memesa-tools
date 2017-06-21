"""
This substituus polco.incore.sh or polco.sh

# TODO: 
- polco memory settings 

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: July 07, 2015
"""

from __future__ import division, print_function, absolute_import

__DEBUG__ = False
input_dir = 'input'
temp_dir = 'temp'
output_dir = 'output'

### Import dependencies ###

import os,shutil,filecmp,sys,re,gzip

import numpy as np

try:
    import sympy
except ImportError:
    print("Error: SymPy is not installed")
    sys.exit()
        
### Check file paths ###

if not os.path.exists(input_dir):
    print("Error: input directory does not exist")
    sys.exit()
   
if not os.path.exists(temp_dir):
    os.mkdir(temp_dir)
   
if not os.path.exists(output_dir):
    os.mkdir(output_dir)    
    
if not os.path.exists('rays'):
    os.mkdir('rays')         

### Set CLASSPATH ###

class TimoTimo():
    def __init__(self,model_name,mode='incore',fixed_identifiers = ["R_module_input","R_module_output"]):
        self.model_name = model_name       
        
        self.fixed_identifiers = fixed_identifiers
        self.model_file = "{0:s}.noinf_r.ine".format(model_name)
        self.model_out_file = "{0:s}.noinf_r.ine.reduced".format(model_name)         
        
        self.getFixedIndices()
        self.removeFixedIndices()
        self.strip2matrix()
        self.runPolco(mode)
        self.countVerticesRays() 
        self.checkScaling()        

        self.extractVerticesRays()
        self.recoverVertices()
        self.checkVertices()
        
        self.recoverRays()
        self.combineVerticesRays()
        self.compressOutput()
        
        
    def getFixedIndices(self):
        """
        We assume that the module input and module output are the only fixed reactions in the system.
        These fixed reactions will be removed before Polco enumerates the vertices and rays. 
        
        alternative approach: provide a file with FVA results
        """
        self.fixed_indices = []            
        file_name = "{0:s}.noinf_r.columns.txt".format(self.model_name)
        if not os.path.isfile(os.path.join(input_dir,file_name)):
           print("Error: File {0:s} not found in directory {1:s}".format(file_name,input_dir))
           sys.exit()
        with open(os.path.join(input_dir,file_name),'r') as infile:
           for line in infile:                         
               if any(ext in line for ext in self.fixed_identifiers): 
                   sline = line.strip().split(',')                    
                   self.fixed_indices.append(sline[0])
         
        self.fixed_indices = map(int,self.fixed_indices)
        self.fixed_indices.sort(reverse=True)
        
    
    def removeFixedIndices(self):
        """ Remove the fixed indices row by row """  
        self.nrows = False
        self.ncols = False
        with open(os.path.join(temp_dir,self.model_out_file), 'w') as outfile:
            with open(os.path.join(input_dir,self.model_file),'r') as infile:
                for line in infile:
                    m = re.match("\d+ +\d+ +rational",line)
                    if m:
                        dimensions = re.findall("\d+",m.group(0))
                        self.nrows = int(dimensions[0])
                        self.ncols = int(dimensions[1])             
                        #print(self.nrows,self.ncols)
                        ncols_reduced = self.ncols - len(self.fixed_indices)
                        outfile.write("H-representation\nbegin\n{0:d} {1:d} rational\n".format(self.nrows,ncols_reduced)) 
                    elif self.nrows and self.ncols:                    
                        sline = line.strip().split(' ')
                        if len(sline) == self.ncols:       
                            sline_sympy = map(sympy.Rational,sline)        
                            for j in self.fixed_indices:                  
                               sline_sympy[0] += sline_sympy[j+1]
                               sline_sympy.pop(j+1)                                                     
                            outfile.write(" ".join(map(str,sline_sympy)))     
                            outfile.write("\n")    
                        else:
                            break               
            outfile.write("end")     
              
    
    def strip2matrix(self):
        """ strip matrix """     
        with open(os.path.join(temp_dir,"temple2"), 'w') as outfile:
            with open(os.path.join(input_dir,self.model_file),'r') as infile:    
                n=1
                for line in infile:
                    sline = line.strip().split(' ')      
                    if len(sline) == self.ncols:                     
                        if n == self.nrows + 1:                             
                            sline[0]='-1'                
                                
                        outfile.write(" ".join(sline))
                        outfile.write("\n") 
                        n+=1
                        

    def runPolco(self,mode='incore',memory = ''): # TODO: memory settings  
        """    
        run polco to enumerate extreme rays of polyhedral cones (http://www.csb.ethz.ch/tools/polco)
        
        Input:
         - *mode* (incore, outcore)
         - *memory* (string)
        """    
        if mode.lower() == 'outcore':
            os.system("java -Xmx3g -Xms3g -Xmn2800M -jar polco.jar -sortinitial LexMin -adj rankup-modpi-outcore -memory out-core -kind cdd -in {0:s} -out text {1:s}.unpol -tmpdir ./rays".format(os.path.join(temp_dir,self.model_out_file),os.path.join(temp_dir,self.model_file))) 

            #os.system("java -Xmx7g -Xms7g -Xmn6800M -jar polco.jar -sortinitial LexMin -adj rankup-modpi-outcore -memory out-core -kind cdd -in {0:s} -out text {1:s}.unpol -tmpdir ./rays".format(os.path.join(temp_dir,self.model_out_file),os.path.join(temp_dir,self.model_file)))     
        else: # assume incore
            os.system("java -Xmx3g -Xms3g -Xmn2800M -jar polco.jar -sortinitial LexMin -kind cdd -in {0:s} -out text {1:s}.unpol -tmpdir ./rays".format(os.path.join(temp_dir,self.model_out_file),os.path.join(temp_dir,self.model_file))) 
       
            #os.system("java -Xmx7g -Xms7g -Xmn6800M -jar polco.jar -sortinitial LexMin -kind cdd -in {0:s} -out text {1:s}.unpol -tmpdir ./rays".format(os.path.join(temp_dir,model_out_file),os.path.join(temp_dir,model_file)))                 
        
        if __DEBUG__: print("Unsorted vertices are now in {0:s}.unpol".format(self.model_file))                    
        
    
    def countVerticesRays(self):
        """ count the number of vertices and rays """         
        polco_delimiter = '\t' 
        if __DEBUG__: print("The vertices will now be in {0:s}.pol and {0:s}.unpol.".format(self.model_file))
        if __DEBUG__: print("Tidying up tabs...")
        self.nrays = 0
        self.nvertices = 0
        with open(os.path.join(temp_dir,"{0:s}.unpol.tidy".format(self.model_file)), 'w') as outfile:
            with open(os.path.join(temp_dir,"{0:s}.unpol".format(self.model_file)),'r') as infile:
                for line in infile:
                    if line.startswith("0"):
                        self.nrays +=1
                    elif line.startswith("1"):
                        self.nvertices +=1

                    fields = line.split(polco_delimiter)
                    outfile.write(' '.join(fields))

        print("---------------------------------------------------------------------------------")
        print("Number of rays:\n{0:d}".format(self.nrays))    
        print("Number of vertices:\n{0:d}".format(self.nvertices))
        print("---------------------------------------------------------------------------------")
        
    
    def checkScaling(self):
        """ check whether Polco has used scaling """              
        if __DEBUG__: print("Going to check whether Polco has used scaling...looking for values that don't start with a 0 or a 1 (looking in {0:s}.unpol.tidy)".format(self.model_file))
        if __DEBUG__:  print("--- begin diff ---")
        filepath_in = os.path.join(temp_dir,"{0:s}.unpol.tidy".format(self.model_file))
        shutil.copy(filepath_in,os.path.join(temp_dir,"danger"))
        IsScaling = False
        with open(filepath_in,'w') as outfile:
            with open(os.path.join(temp_dir,"danger"),'r') as infile:
                for line in infile:
                   sline = line.strip().split(' ')               
                   if sline[0] in ['0','1']: # OK
                       outline = sline
                   else:                     # Scaling: divide by first entree
                       IsScaling = True
                       sline_sympy = np.array(map(sympy.Rational,sline))
                       outline = map(str,sline_sympy/sline_sympy[0])
                   
                   outfile.write(" ".join(outline))
                   outfile.write("\n")   
                        
        #os.system("java PolcoScale danger > {0:s}".format(filepath_in))         
        if __DEBUG__:  
            if IsScaling:
                print("Polco used scaling ... ")
            else:
                print("Polco used no scaling ...")    
            print("--- end diff ---")        
                   
    
    def extractVerticesRays(self):
        """ put vertices and rays in different files """              
        if __DEBUG__:  print("Going to extract the vertices and rays into different files ...")        
        outfile_ver = open(os.path.join(temp_dir,"{0:s}.unpol.tidy.onlyver".format(self.model_file)), 'w')
        outfile_ray = open(os.path.join(temp_dir,"{0:s}.unpol.tidy.onlyray".format(self.model_file)), 'w')
        with open(os.path.join(temp_dir,"{0:s}.unpol.tidy".format(self.model_file)), 'r') as infile:
            for line in infile:
                if not line.startswith("0"):
                    fields = line.split(' ')
                    if '' in fields:
                        fields.remove('')
                    outfile_ver.write(' '.join(fields[1:]))
                else: 
                    fields = line.split(' ')
                    outfile_ray.write(' '.join(fields[1:]))
        outfile_ver.close()
        outfile_ray.close()
        
    
    def recoverVertices(self):
        """ recover vertices of the original system (i.e. put fixed fluxes back in)  """
        if __DEBUG__: print("Going to recover the vertices of the original system (i.e. put fixed fluxes back in)")        
        
        self.fixed_indices.sort()        
        with open(os.path.join(temp_dir,"{0:s}.unpol.tidy.restored.onlyver".format(self.model_file)), 'w') as outfile:
            with open(os.path.join(temp_dir,"{0:s}.unpol.tidy.onlyver".format(self.model_file)), 'r') as infile:       
                for line in infile:     
                    fields = line.strip().split(' ')                                       
                    for j in self.fixed_indices:
                        fields.insert(j,'1') # TODO: this only works for the input and output reaction. For others, we need FVA data
                    outfile.write(' '.join(fields))
                    outfile.write('\n')
                    
                    
    def checkVertices(self,ncheck=10):
        """
        Check if the enumerated vertices are vertices in the original space by determining:
          - if they satisfy the constraints (i.e. are inside the polytope)
          - if the set of the tight inequalities has a full rank
        """
        if __DEBUG__: print("Finished converting the points. Am now going to check that the first 10 !VERTICES! (not the rays) are all vertices in the original space...")
        if __DEBUG__: print("Any points that are not vertices (or not even in the polytope) will elicit the message {F}AKE VERTEX!")
       
        #MIDROWS = self.nrows + 1
        #MIDCOLS = self.ncols
        MIDFILE = os.path.join(temp_dir,"{0:s}.unpol.tidy.restored.onlyver".format(self.model_file))
        #os.system("java CheckVertices {0:s} {1:d} {2:d} {3:s} | tee {4:s}.report.ver".format(os.path.join(temp_dir,"temple2"),MIDROWS,MIDCOLS,MIDFILE,os.path.join(temp_dir,self.model_file)))
        
        print("***************************************")
        Arr_b = []
        Arr_A = []
        Arr_Ab = []                                
  
        ### Import temple2 which contains the Ax >= b data 
        with open(os.path.join(temp_dir,"{0:s}.report.ver".format(self.model_file)),"w") as outfile:
            with open(os.path.join(temp_dir,"temple2"),'r') as infile:
                for line in infile:
                   sline = line.strip().split(" ")
                   sline_sympy = map(sympy.Rational,sline)        
                   Arr_b.append(self.Negate(sline_sympy[0]))
                   Arr_A.append(sline_sympy[1:])
                   Arr_Ab.append(sline_sympy)

            Arr_b = np.array(Arr_b)
            Arr_A = np.array(Arr_A)
            Arr_Ab = np.array(Arr_Ab)

            rows = Arr_A.shape[0]
            cols = Arr_A.shape[1]

            #### Import vertices
            with open(MIDFILE,'r') as infile:       
                for n,vertex in enumerate(infile):        # for each vertex in infile    
                   if n==ncheck:
                      break
                   svertex = vertex.strip().split(" ")
                   svertex_sympy = map(sympy.Rational,svertex)    
                   Arr_v = np.transpose(svertex_sympy) 
                   print("* Computing Av")
                   outfile.write("* Computing Av\n")
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
                          outfile.write("* IMPORTANT: This 'vertex' was not even in the polytope!\n* This is a FAKE VERTEX\n")            
                          fake_vertex=True        
                          break
                      i+=1

                   if not fake_vertex:  
                       ### Determine the rank of the submatrix which must be full rank
                       print("* So {0:d} of the inequalities were tight. Checking the rank of this set ...".format(ntight))
                       outfile.write("* So {0:d} of the inequalities were tight. Checking the rank of this set ...\n".format(ntight))
                       #bgoli 2017 here I do an explicit cast to float64 for the svd evaluation
                       Arr_sub_Ab = np.array(Arr_sub_Ab)
                       u_, s_, v_ = np.linalg.svd(Arr_sub_Ab.astype('float64'))
                       matrix_rank = np.sum(s_ > 1e-10)
                       print("* The submatrix has rank {0:d}".format(matrix_rank))
                       outfile.write("* The submatrix has rank {0:d}\n".format(matrix_rank))                      
                       if matrix_rank == cols:  # full rank
                           print("* This is a REAL VERTEX")           
                           outfile.write("* This is a REAL VERTEX\n")           
                       else: 
                           print("* This is a FAKE VERTEX") 
                           outfile.write("This is a FAKE VERTEX\n")           
        

        if __DEBUG__: print("Am putting the 'vertices' to the file {0:s}.ver, a copy of the vertex-checking is in {0:s}.report.ver".format(self.model_file))

        filepath_in = os.path.join(temp_dir,"{0:s}.unpol.tidy.restored.onlyver".format(self.model_file))
        filepath_out = os.path.join(temp_dir,"{0:s}.ver").format(self.model_file)
        shutil.copy(filepath_in,filepath_out)


    def recoverRays(self):
        """ recover rays of the original system (i.e. put fixed fluxes back in) """
        if __DEBUG__: print("Going to convert the rays (and linealities) in the column stripped system, back into the original system")                    
        
        self.fixed_indices.sort()
        with open(os.path.join(temp_dir,"{0:s}.unpol.tidy.onlyray".format(self.model_file)), 'r') as infile:
           with open(os.path.join(temp_dir,"{0:s}.unpol.tidy.restored.onlyray".format(self.model_file)), 'w') as outfile:
                for line in infile:     
                    fields = line.split(' ')                            
                    for j in self.fixed_indices:
                        fields.insert(j,'0')         
                    outfile.write(' '.join(fields))   
       
        if __DEBUG__: print("Am putting the 'rays' to the file {0:s}.rays.".format(self.model_file))
        shutil.copy(os.path.join(temp_dir,"{0:s}.unpol.tidy.restored.onlyray".format(self.model_file)),os.path.join(temp_dir,"{0:s}.ray".format(self.model_file)))
        

    def combineVerticesRays(self):
        if __DEBUG__: print("Will combine (lineality space) rays and vertices into the file{0:s}.all".format(self.model_file))

        filepath_out = os.path.join(temp_dir,"{0:s}.all".format(self.model_file))
        with open(filepath_out, 'w') as outfile:
             with open(os.path.join(temp_dir,"{0:s}.ray".format(self.model_file)), 'r') as infile:
                outfile.write("* Rays ({0:d} vectors)\n".format(self.nrays))
                for line in infile:                          
                    outfile.write(line)    
             with open(os.path.join(temp_dir,"{0:s}.ver".format(self.model_file)), 'r') as infile:
                outfile.write("* Vertices ({0:d} vectors)\n".format(self.nvertices))
                for line in infile:                             
                    outfile.write(line)                      


    def compressOutput(self):        
        f = gzip.open(os.path.join(output_dir,"{0:s}.all.gz".format(self.model_file)),'wb')
        with open(os.path.join(temp_dir,"{0:s}.ray".format(self.model_file)), 'r') as infile:
            f.write("* Rays ({0:d} vectors)\n".format(self.nrays))
            for line in infile:                          
                f.write(line)    
        with open(os.path.join(temp_dir,"{0:s}.ver".format(self.model_file)), 'r') as infile:
            f.write("* Vertices ({0:d} vectors)\n".format(self.nvertices))
            for line in infile:                             
                f.write(line)            
        f.close()        
        
    def Negate(self,x):
        if x:
            x = -x
        return(x)             

TimoTimo("toy_model.1_split")
TimoTimo("toy_model.2_split")
TimoTimo("toy_model.3_split")
    
    

    
    
