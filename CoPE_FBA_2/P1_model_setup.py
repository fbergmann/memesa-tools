"""
Get Model Information
=====================

This is the first script of the CoPE-FBA 2.0 pipeline. I will ask for the *model_name* which must be in the same directory as this script. You can also specify the *model_name* as first and only argument.

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

import os,shutil,sys

try:
    import libsbml
except ImportError:
    print('\nI require libSBML with Python bindings: http://sbml.org/Software/libSBML\n')    
cDir = os.getcwd()

if len(sys.argv) == 2:  
   print('Number of arguments:', len(sys.argv), 'arguments.')
   print( 'Argument List:', str(sys.argv)) 
   model_name = sys.argv[1]
else:
    try: input = raw_input
    except NameError: pass
    model_name = input('Specify the model filename please:\n')
    
if model_name.endswith('.xml'):
    model_name = model_name.replace('.xml','')
       
if not os.path.exists('data'): os.mkdir('data')
if not os.path.exists(os.path.join('data',model_name)): os.mkdir(os.path.join('data',model_name))
if not os.path.exists(os.path.join('data',model_name,'flux_modules')): os.mkdir(os.path.join('data',model_name,'flux_modules'))
if not os.path.exists(os.path.join('data',model_name,'models')): os.mkdir(os.path.join('data',model_name,'models'))
if not os.path.exists(os.path.join('data',model_name,'models','sbml')): os.mkdir(os.path.join('data',model_name,'models','sbml'))
if not os.path.exists(os.path.join('data',model_name,'cope_fba')): 
    os.mkdir(os.path.join('data',model_name,'cope_fba'))
    os.mkdir(os.path.join('data',model_name,'cope_fba','subnetworks'))
    os.mkdir(os.path.join('data',model_name,'cope_fba','subnetworks','vertex'))
    
print("Info: directories generated")

D = libsbml.readSBML("{0:s}.xml".format(model_name))
M = D.getModel()
if D.getLevel() == 3:
    SBML = 3
    shutil.copy2("{0:s}.xml".format(model_name), os.path.join('data',model_name,'models','sbml'))
    print("Info: your model is written in SBML L3")
elif "<fba:listOfConstraints>" in M.getAnnotationString():
    SBML = 2 
    shutil.copy2("{0:s}.xml".format(model_name), os.path.join('data',model_name,'models','sbml'))
    print("Info: your model is written in SBML L2")
else:
    print('I suspect this is a COBRA file, attempting to convert')
    try:
        import pyscescbm as cbm
    except ImportError:
        print('\nI require PySCeS-CBMPy: http://cbmpy.sourceforge.net\n')
    cobramod = cbm.readCOBRASBML("{0:s}.xml".format(model_name))
    cbm.writeSBML3FBC(cobramod, os.path.join('data',model_name,'models','sbml','{0:s}.xml'.format(model_name)))
    SBML = 3
    print("Info: we converted your model in SBML L3")
    
try:
    import numpy
except ImportError:
    print('\nI require NumPy: http://www.numpy.org\n')    

try:
    import sympy    
except ImportError:
    print('\nI require SymPy: http://sympy.org/en/index.html\n')

try:
    import h5py
except ImportError:
    print('\nI require H5Py http://www.h5py.org/\n')  
