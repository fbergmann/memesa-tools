: '
Run Analysis Example script
==============================

In this example, we illustrate how you can quickly analyze the enumerated vertices

Make sure that the corresponding SBML file is in the home folder of CoPE-FBA 2.0

Before we can execute the sh script, we change the access permissions: chmod 777 run_analysis_example.sh

The analysis phase consists of scripts P8-P11 and we can pass the following arguments:
python P8_vertex_translate.py -m *model_name* -l *sbml_level*
python P9_vertex_statistics.py -m *model_name* -l *sbml_level*
python P10_secondary_optimization.py -m *model_name* -l *sbml_level*
python P11_ray_statistics.py -m *model_name* -l *sbml_level*

Some scripts (P9 and P10) accept additional arguments. Please read the documentation of these scripts carefully.

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, B. Teusink, and F.J. Bruggeman.
'

# 1. toy model
python P8_vertex_translate.py -m toy_model -l 3
python P9_vertex_statistics.py -m toy_model -l 3
python P10_secondary_optimization.py -m toy_model -l 3
python P11_ray_statistics.py -m toy_model -l 3

# 2. iTM686.light
python P8_vertex_translate.py -m iTM686.light -l 2
python P9_vertex_statistics.py -m iTM686.light -l 2
python P10_secondary_optimization.py -m iTM686.light -l 2
python P11_ray_statistics.py -m iTM686.light -l 2

# 3. Ecoli_iAF1260_noox.glc
python P8_vertex_translate.py -m Ecoli_iAF1260_noox.glc -l 2
python P9_vertex_statistics.py -m Ecoli_iAF1260_noox.glc -l 2
python P10_secondary_optimization.py -m Ecoli_iAF1260_noox.glc -l 2
python P11_ray_statistics.py -m Ecoli_iAF1260_noox.glc -l 2
