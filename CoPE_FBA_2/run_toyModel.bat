@echo off
Setlocal EnableDelayedExpansion
SET BASE_DIR=%~dp0

SET _PATH=%PATH%
SET _CP=%CLASSPATH%
SET PYTHON_DIR=C:\Anaconda2_32;C:\Anaconda2_32\script
SET BASH_DIR=C:\Program Files\Git\bin;C:\Program Files\Git\usr\bin
SET JAVA_DIR=C:\Program Files\Java\jdk1.8.0_131\bin

SET PATH=%PYTHON_DIR%;%BASH_DIR%;%JAVA_DIR%;%PATH%
SET CLASSPATH=%BASE_DIR%/rational_fba;%BASE_DIR%/rational_fba/jlinalg_0.5.jar;%CLASSPATH%

echo P1
python P1_model_setup.py toy_model

echo P2
python P2_find_duplicate_reactions.py -m toy_model -b 1000 -l 3

echo P3
python P3_get_F_modules.py -m toy_model -l 3

echo P4
python P4_SBML2H-format.py -m toy_model -b 1000 -l 3

echo P5
sh ./P5_do_rational_FBA.sh toy_model

echo P6
python P6_stripModel2Modules.py -m toy_model -b 1000 -l 3


:: restore path 
SET PATH=%_PATH%
SET CLASSPATH=%_CP%