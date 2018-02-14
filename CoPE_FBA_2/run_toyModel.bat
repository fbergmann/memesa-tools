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
SET MODEL=toy_model
SET INFINITY=1000
SET LEVEL=3

echo P1
python P1_model_setup.py %MODEL%

echo P2
python P2_find_duplicate_reactions.py -m %MODEL% -b %INFINITY% -l %LEVEL%

echo P3
python P3_get_F_modules.py -m %MODEL% -l %LEVEL%

echo P4
python P4_SBML2H-format.py -m %MODEL% -b %INFINITY% -l %LEVEL%

echo P5
cd rational_fba
sh ./rational_FBAw.sh %MODEL%
cd ..

echo P6
python P6_stripModel2Modules.py -m %MODEL% -b %INFINITY% -l %LEVEL%


:: restore path 
SET PATH=%_PATH%
SET CLASSPATH=%_CP%