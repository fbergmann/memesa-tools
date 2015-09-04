CoPE-FBA 2.0: Fast Enumeration of the optimal solution space for genome-scale metabolic models 
==============================================================================================

CoPE-FBA (Comprehensive Polyhedron Enumeration Flux Balance Analaysis) characterizes the optimal solution space in terms of vertices, linealities and rays. We present CoPE-FBA 2.0, an updated pipeline which allows for fast enumeration of the set of vertices and rays. We split reversible reactions, so linealities do not exist. The set of vertices now represents all non-decomposable flux vectors in the optimal flux space.

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, F.J. Bruggeman, and B. Teusink

Pipeline usage
--------------

We present a collection of scripts that forms the pipeline to enumerate the set of vertices and rays. This pipeline exists of of three parts:

(1) file preparation (P1-P6)
(2) data generation (P7, Polco)
(3) data analysis (P8-P11)

To make this pipeline more user-friendly, we developed several shell script files which perform the three parts of the pipeline.  Alternatively, one can run each file of the pipeline after each other. Both the file preparation and the data generation steps are essential, but users are free to perform any kind of analysis. We recommend using the pipeline up to at least P8. 

The arguments (hereafter USER_INPUT) can be set in the scripts or parsed as arguments. See the documentation in the scripts for more information.

File description
----------------

(P1) Setup file directories with P1_model_setup.py
USER_INPUT: *model_name* 

Protein costs can be incorporated. Provide a csv file named according to the *model_name* and the type of cost used. Put this csv file in the protein_costs subdirectory. Example: toy_model/protein_costs/toy_model.avg.costs.csv, toy_model/protein_costs/toy_model.max.costs.csv

(P2) Scan for duplicate reactions with P2_find_duplicate_reactions.py
USER_INPUT: *model_name*, *inf_bound*, and *sbml_level*.
Remove all modeling duplicates, artefacts and isozymes. Biologically relevant ones, e.g. atp maintenance vs nucleotide biosynthesis, should not be removed!

(P3) Enumerate the F-modules with P3_get_F_modules.py.
(Essential) USER INPUT: *model_name* and *sbml_level*

(P4) Generate the H-representation (.ine extension) with P4_SBML2H-format.py.
USER_INPUT: *model_name*, *inf_bound*, and *sbml_level*.

(P5) Perform rational FBA from the H-representation generated at the previous step
USER_INPUT: *model_name* (user input must be provided as an argument --> ./P5_rational_FBA.sh *model_name*)

(P6) Use P6_stripModel2Modules.py to generate subnetwork models with rational input-output (and a dummy variable) in different formats: H-representation (.ine), SBML (.xml), MATLAB (.m).
(Essential) USER INPUT: *model_name*, *inf_bound*, and *sbml_level*

(P7) Enumerate vertices and rays with Polco 
USER_INPUT: *H-representation* (the .ine files generated at step 6.) 

We recommend users of this pipeline to store the polco directory somewhere in the home folder.
       
Polco generates an archive file format which should be extracted in the following directory:
/*model_name*/cope_fba/subnetworks/vertex       
  
(P8) Translate the enumerated sets of vertices and rays with P8_vertex_translate.py 
USER INPUT: *model_name* and *sbml_level*

(P9) Determine whole model vertex statistics with P9_vertices_statistics.py
(Essential) USER INPUT: *model_name* and *sbml_level*
Use *DO_FULL_ENUMERATE*=TRUE to get the vertex length, cost etc. for every vertex. Warning, this can take a significant amount of time if there are many vertices

(P10) Use MILP to double check the outcomes of the secondary objectives with P10_secondary_optimization.py
(Essential) USER INPUT: *model_name* and *sbml_level*

(P11) Get information about the rays in the F-modules with 11_ray_statistics.py 
(Essential) USER INPUT: *model_name* and *sbml_level*
