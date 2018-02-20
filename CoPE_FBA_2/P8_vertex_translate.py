"""
Translate Polco output to HDF5 and CSV files
============================================

User input: 
- *model_name* (string)
- *sbml_level* (int)

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

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 25, 2015
"""

from __future__ import division, print_function, absolute_import

import os, time,sys,getopt
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
ALTLOADFILE = None

# Set userdata here #

model_name = 'toy_model'
sbml_level = 3

# optional user data #

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

# ALTLOAD do not do .ine to hdf5 conversion simply load the specified hdf5 file
ALTLOAD = False
# use the *.columns.txt file from the Steven output to cross check flux order
USE_COLUMN_CROSSCHECK = True
# rounding off digits
rFact = 15
# write vertex_variable_flux.csv
FIND_VARIABLE_FLUXES_IN_VERTICES = True # needed for module analysis
#  write vertex_array.csv
WRITE_VERTEX_ARRAY = False # not needed
WRITE_VERTEX_ARRAY_FULL = False # once only small models
# memory efficient algorithms
BIG_FILE = True
# end userdata here #

import numpy, sympy, h5py, gc
gc.enable()

# bgoli 2017-06-21
# CBMPy compatability changes
import cbmpy
__DEBUG__ = cbmpy.CBConfig.__CBCONFIG__['DEBUG']
__version__ = cbmpy.CBConfig.current_version()
import cbmpy.fluxmodules.fluxmodules as fmod
from cbmpy import CBRead, CBWrite, CBTools, CBMultiCore
from cbmpy import CBSolver as slv

model_dir = os.path.join(cDir, 'data', model_name,'models_subnetwork','sbml')
H_format_dir = os.path.join(cDir,'data', model_name,'models_subnetwork','h-format')
work_dir = os.path.join(cDir, 'data',model_name,'cope_fba','subnetworks', 'vertex')

def writeVertexArray(rng, model_name, mod, vertList, vertex_arr_status, vertex_arr_min, vertex_arr_max):
    fname = os.path.join(work_dir,'{0:s}.vertex_array.csv'.format(model_name))   
    FO = open(fname, 'w')
    cntr = 1
    fcntr = 0
    col_l = len(vertex_arr_status)
    row_l = vertList.shape[0]
    print( 'vl', row_l, col_l, vertList.shape[1], rng)
    Btime = time.time()
    for col in range(*rng):
        if vertex_arr_status[col] == 0.0:
            print(time.strftime('%d-%H:%M'), '%3.2f' % (float(cntr)/float(realFVA_vari_len)*100), '%2.2f' % (float(time.time()-Btime)/60.0), mod.N.col[col], col)
            cntr += 1
            R = mod.getReaction(mod.N.col[col])
            if 'SUBSYSTEM'  in R.annotation:
                subsystem = '"{0:s}"'.format(R.annotation['SUBSYSTEM'])
            else:
                subsystem = ''
            if 'Equation' in R.annotation:
                equation = '"{0:s}"'.format(R.annotation['Equation'])
            else:
                equation = ''
            FO.write('"{0:s}","{1:s}",{2:s},{3:s},{4},{5},'.format(R.getPid(), R.name, equation, subsystem, vertex_arr_min[col],vertex_arr_max[col]))
            outStr = ''
            rcntr = 0
            for row in range(row_l):
                rcntr += 1
                if row != row_l-1:
                    if vertList[row][col] == 0.0:
                        outStr += '0.0,'
                    else:
                        outStr += '{0},'.format(round(vertList[row][col],rFact))
                else:
                    if vertList[row][col] == 0.0:
                        outStr += '0.0\n'
                    else:
                        outStr += '{0}\n'.format(round(vertList[row][col],rFact))
                if rcntr == 5000:
                    FO.write(outStr)
                    rcntr = 0
                    outStr = ''
            FO.write(outStr)
            FO.flush()
            outStr = ''
    FO.close()    
    return fname

if ALTLOAD:
    ALTLOADFILE = os.path.join(work_dir, ALTLOADFILE)
    assert os.path.exists(ALTLOADFILE)    

if not os.path.exists(work_dir):
    os.mkdir(work_dir)   

if os.listdir(work_dir) == []:
    raise UserWarning("No files to analyze in directory {0:s}".format(work_dir))

for file_in in os.listdir(work_dir):       
    if file_in.endswith('noinf_r.ine.all'):
        subnetwork_name = file_in.replace('.noinf_r.ine.all','') 
        model_file = "{0:s}.xml".format(subnetwork_name)
        vertex_file = os.path.join(work_dir, '{0:s}.noinf_r.ine.all'.format(subnetwork_name) )
        rfva_file = os.path.join(work_dir, '{0:s}.noinf_r.ine.opt.fva'.format(subnetwork_name) )   
        if sbml_level == 3:
            try:
                cmod = CBRead.readSBML3FBC(model_file, work_dir=model_dir)  
            except:        
                cmod = CBRead.readSBML2FBA(model_file, work_dir=model_dir)      
        else:    
            cmod = CBRead.readSBML2FBA(model_file, work_dir=model_dir)        
        cmod.id = subnetwork_name
        CBTools.addStoichToFBAModel(cmod)
        CBTools.processBiGGchemFormula(cmod)
        
        if USE_COLUMN_CROSSCHECK:
            mismatch = []                        
            F = open(os.path.join(H_format_dir, '{0:s}.noinf_r.columns.txt'.format(subnetwork_name) ), 'r')
            indx_lst = []
            for L in F:
                l = [i.strip() for i in L.split(',')]
                if len(l) == 2:
                    indx_lst.append(l[1])
            F.close()
            for r in range(len(indx_lst)):
                match = bool(indx_lst[r] == cmod.N.col[r])
                print('{0} == {1}: {2}'.format(indx_lst[r], cmod.N.col[r], match))
                if not match: mismatch.append((indx_lst[r], cmod.N.col[r]))
            if len(mismatch) > 0:
                print('\nColumn mismatch ({0:d})\n'.format(len(mismatch)) )
                ##  print mismatch
                os.sys.exit(1)
            else:
                print('\nColumn check successful, have a nice day :-)')
                time.sleep(2)        
                
        TIME_START = time.time()
        # get Real FVA results from Steven
        realFVAdata = CBRead.readSK_FVA(rfva_file)

        for f in range(len(realFVAdata)):
            realFVAdata[f].update({'id' : cmod.N.col[f]})

        realNamesVar = [j['id'] for j in realFVAdata if j['status']=='VARIABLE']
        realFVA_vari_len = len(realNamesVar)
        realNamesFix = [j['id'] for j in realFVAdata if j['status']=='FIXED']
        realFVA_fixed_len = len(realNamesFix)
        print("\nReal FVA reports {0:d} fixed and {1:d} variable fluxes.".format(realFVA_fixed_len, realFVA_vari_len))        
        
        FLUX_NAMES = tuple(cmod.N.col)

        # parse vertex file to array and symb_array
        if not ALTLOAD:
            HD5datF = CBRead.readSK_vertex(vertex_file, bigfile=True, fast_rational=True, nformat='%.15f', compression=None, hdf5file=os.path.join(work_dir, subnetwork_name))
        else:
            HD5datF = ALTLOADFILE
        HD5dat = h5py.File(HD5datF,'r')

        print('\nUsing HDF5 datafile: {0:s}\n'.format(HD5datF) )

        HAVE_RAYS = True
        HAVE_LIN = True

        Tlist = []
        try:
            for l in range(HD5dat['data/rays'].shape[0]):
                Tlist.append(HD5dat['data/rays'][l])
            rayList = Tlist
        except:
            print('No Rays')
            rayList = []
            HAVE_RAYS = False
        Tlist = []
        try:
            for l in range(HD5dat['data/lin'].shape[0]):
                Tlist.append(HD5dat['data/lin'][l])
            linBasis = Tlist
        except:
            print('No Lineality vectors')
            linBasis = []
            HAVE_LIN = False
        del Tlist
        print('Rays & Linealities {0:d} {1:d}'.format( len(rayList), len(linBasis) ) )

        rayF = open(os.path.join(work_dir,'{0:s}.rays.csv'.format(subnetwork_name)),'w')
        if HAVE_RAYS:
            # find fluxes active in rays
            ray_arr = numpy.array(rayList)
            RAY_FLUXES = []
            del rayList
            #print(ray_arr.shape)
            for c in range(ray_arr.shape[1]):
                if (ray_arr[:,c] > 0.0).any():
                    #print FLUX_NAMES[c], (ray_arr[:,c] > 0.0).any()
                    RAY_FLUXES.append(FLUX_NAMES[c])
            print(RAY_FLUXES)

            # write the rays
            for c in range(len(realFVAdata)):
                if c != len(realFVAdata)-1:
                    rayF.write('{0},'.format(c+1))
                else:
                    rayF.write('{0}\n'.format(c+1))
            for c in range(len(realFVAdata)):
                id = realFVAdata[c]['id']
                if c != len(realFVAdata)-1:
                    rayF.write('{0},'.format(id) )
                else:
                    rayF.write('{0}\n'.format(id) )
            ray_arr.shape
            for r in range(ray_arr.shape[0]):
                for c in range(ray_arr.shape[1]):
                    if c != ray_arr.shape[1]-1:
                        rayF.write('{0},'.format(ray_arr[r,c]) )
                    else:
                        rayF.write('{0}\n'.format(ray_arr[r,c]) )
            del ray_arr
        rayF.flush()
        rayF.close()

        linbF = open(os.path.join(work_dir,'{0:s}.lineality_basis.csv'.format(subnetwork_name) ),'w')
        if HAVE_LIN:
            # find fluxes active in lineality space
            linbasis_arr = numpy.array(linBasis)
            del linBasis
            #print(linbasis_arr.shape)

            LIN_FLUXES = []
            for c in range(linbasis_arr.shape[1]):
                if (linbasis_arr[:,c] > 0.0).any():
                    #print FLUX_NAMES[c], (linbasis_arr[:,c] > 0.0).any()
                    LIN_FLUXES.append(FLUX_NAMES[c])

            #print(len(LIN_FLUXES))
            print(LIN_FLUXES)

            # write the basis for the lineality space
            for c in range(len(realFVAdata)):
                if c != len(realFVAdata)-1:
                    linbF.write('{0},'.format(c+1) )
                else:
                    linbF.write('{0}\n'.format(c+1) )
            for c in range(len(realFVAdata)):
                id = realFVAdata[c]['id']
                if c != len(realFVAdata)-1:
                    linbF.write('{0},'.format(id) )
                else:
                    linbF.write('{0}\n'.format(id) )
            for r in range(linbasis_arr.shape[0]):
                for c in range(linbasis_arr.shape[1]):
                    if c != linbasis_arr.shape[1]-1:
                        linbF.write('{0},'.format(linbasis_arr[r,c]) )
                    else:
                        linbF.write('{0}\n'.format(linbasis_arr[r,c]) )
            del linbasis_arr
        linbF.flush()
        linbF.close()       
        
        
        # get the vertices as a list/array
        try:
           vertList = HD5dat['data/vertices']
           HAVE_VERTICES = True
        except:
           print('No Vertices')
           HAVE_VERTICES = False
        #print len(vertList)

        if  HAVE_VERTICES:
            # populate matrices with fixed data
            vertex_arr_status = []
            vertex_arr_min = []
            vertex_arr_max = []
            for j in realFVAdata:
                j_min = sympy.Rational('{0}'.format(j['min']) )
                j_max = sympy.Rational('{0}'.format(j['max']) )
                vStat = None
                if j['status']=='FIXED':
                    vStat = 1
                else:
                    vStat = 0
                vertex_arr_status.append(vStat)
                vertex_arr_min.append(j_min.evalf())
                vertex_arr_max.append(j_max.evalf())
            del j_min, j_max

            """
            create a map index file for vertex_array with the following columns:

             - index
             - id
             - name
             - subsystem
             - fixed (1/0)
             - FVA min
             - FVA max

            The last three columns are included in vertex_array as the first 3 rows
            """
            # write vertex_array_header.csv
            mapF = open(os.path.join(work_dir,'{0:s}.vertex_array_header.csv'.format(subnetwork_name) ),'w')
            mapF.write('Index,Id,Name,Subsystem,Fixed,FVAmin,FVAmax,Equation\n')
            for r in range(len(realFVAdata)):
                R = cmod.getReaction(realFVAdata[r]['id'])
                name = '"{0:s}"'.format(R.name)
                if 'SUBSYSTEM' in R.annotation:
                    subsystem = '"{0:s}"'.format(R.annotation['SUBSYSTEM'])
                else:
                    subsystem = ''
                if 'Equation' in R.annotation:
                    equation = '"{0:s}"'.format(R.annotation['Equation'])
                else:
                    equation = ''
                mapF.write('{0:d},{1:s},{2:s},{3:s},{4},{5},{6},{7}\n'.format(r+1, realFVAdata[r]['id'], name, subsystem, vertex_arr_status[r], vertex_arr_min[r], vertex_arr_max[r], equation))
            mapF.close()
            gc.collect()

            #  write vertex_array_full.csv
            if WRITE_VERTEX_ARRAY_FULL:
                CBTools.exportLabelledLinkedList([vertex_arr_status], os.path.join(work_dir,'{0:s}.vertex_array_full.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=False)
                CBTools.exportLabelledLinkedList([vertex_arr_min], os.path.join(work_dir,'{0:s}.vertex_array_full.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=True)
                CBTools.exportLabelledLinkedList([vertex_arr_max], os.path.join(work_dir,'{0:s}.vertex_array_full.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=True)

                # write out chunks
                cntr = 0
                cntr2 = 2000
                while cntr < len(vertList)-50 and cntr2 < len(vertList)-1:
                    CBTools.exportLabelledLinkedList(vertList[cntr:cntr2], os.path.join(work_dir,'{0:s}.vertex_array_full.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=True)
                    cntr = cntr2
                    cntr2 += 2000

                # write out the remaining bit
                for x in range(cntr, len(vertList)):
                    CBTools.exportLabelledLinkedList([vertList[x]], os.path.join(work_dir,'{0:s}.vertex_array_full.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=True)            
                            
            #  write vertex_array.csv
            if WRITE_VERTEX_ARRAY:
                fname = writeVertexArray((0,len(vertex_arr_status)), subnetwork_name, mod, vertList, vertex_arr_status, vertex_arr_min, vertex_arr_max)
            # find constant fluxes in vertices export as (write vertex_variable_fluxes.csv)
            if FIND_VARIABLE_FLUXES_IN_VERTICES:
                vrl = len(vertList)
                vcl = len(vertList[0])
                ZERO_CUT = 1.0e-10
                #print vrl, vcl
                BoolRes = numpy.zeros(vcl,'bool')
                INITIAL_ROW = vertList[0].copy()

                cntr = 0
                for r in range(vrl):
                    BoolDiff = (numpy.absolute(vertList[r] - INITIAL_ROW) >= ZERO_CUT)
                    for B in range(len(BoolDiff)):
                        if BoolDiff[B]:
                            BoolRes[B] = True
                    if cntr == 1000:
                        print('Processing vertex {0:d} of {1:d}'.format(r, vrl))
                        cntr = 0
                    else:
                        cntr += 1

                VarFluxes = []
                for f in range(len(BoolRes)):
                    if BoolRes[f]:
                        VarFluxes.append(FLUX_NAMES[f])

                print('Variable fluxes:', len(VarFluxes))
                print(VarFluxes)
                CBTools.exportLabelledLinkedList([VarFluxes], os.path.join(work_dir,'{0:s}.vertex_variable_fluxes.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=False) # EvP changed format into fmt
                CBTools.exportLabelledLinkedList([list(FLUX_NAMES), BoolRes.tolist()], os.path.join(work_dir,'{0:s}.vertex_variable_fluxes_all.csv'.format(subnetwork_name) ), sep=',', fmt='%s', appendlist=False)# EvP changed format into fmt


        TIME_END = time.time()        
        print('\n\nTime taken to analyse {0:s} data: {1:2.2f} min\n'.format(subnetwork_name, (TIME_END-TIME_START)/60.0))
