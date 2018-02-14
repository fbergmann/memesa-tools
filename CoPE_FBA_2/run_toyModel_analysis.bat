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



echo P8
python P8_vertex_translate.py -m %MODEL% -l %LEVEL%
echo P9
python P9_vertex_statistics.py -m %MODEL% -l %LEVEL%
echo P10
python P10_secondary_optimization.py -m %MODEL% -l %LEVEL%

echo P11
python P11_ray_statistics.py -m %MODEL% -l %LEVEL%


:: restore path 
SET PATH=%_PATH%
SET CLASSPATH=%_CP%