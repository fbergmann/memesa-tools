## Fork from: SystemsBioinformatics/memesa-tools
This fork makes it possible to run the pipeline on Windows. In order to make this work, I recompiled the esolver using cygwin by checking out: 

* <https://github.com/jonls/qsopt-ex>

Other than that I've created `run_polco_from_name.py`, to make it easier running step P7 of the pipeline. This script just takes the arguments: 

	model-name <input_dir> <output_dir> <filter>

where the `model-name` is mandatory and the other arguments optional, but if given they have to occur in the correct order. the `input_dir` is usually: 

	../data/%MODEL%/models_subnetwork/h-format

the `output_dir`: 

	../data/%MODEL%/cope_fba/subnetworks/vertex

the `filter` argument makes it possible to restrict as to what files ought to be found. the filter string is directly passed into `glob` that matches all files: `*{filter}_split.noinf`. 

Example executions of this process are in this repository for invoking the toy model. The steps are: 

* preparation p1-p6: `run_toyModel.bat`
* enumeration p7: `run_toyModel_p7.bat`
* analysis P8-P11: `run_toyModel_analysis.bat`

### Dependencies
To make this work on Windows the following dependencies need to be installed: 

* Anaconda, including cplex and cbmpy
* Git for Windows, (unix subsystem, shell)
* Java

simply modify the variables `PYTHON_DIR`, `BASH_DIR` and `JAVA_DIR` on the top of the batch scripts.  

### Apply for new models
* copy the model into the directory, 
* copy the toy model batch files, renaming them to match your model
* edit those batch files, replacing the variables `MODEL`, `INFINITY` and `LEVEL`. 

---
2/14/2018 3:21:23 PM  Frank T. Bergmann