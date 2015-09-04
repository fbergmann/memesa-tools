: '
Run Preparation Example script
==============================

In this example, we illustrate how you can quickly generate the H-representation files (.ine) which are used by QSopt_EX (exact LP) and Polco (vertex enumeration).

Make sure that the corresponding SBML file is in the home folder of CoPE-FBA 2.0

Before we can execute the sh script, we change the access permissions: chmod 777 run_preparation_example.sh

The preparation phase consists of scripts P1-P6 and we can pass the following arguments:
python P1_model_setup.py *model_name*
python P3_SBML2H-format.py -m *model_name* -b *inf_bound* -l *sbml_level*
sh ./P4_do_rational_FBA.sh *model_name*
python P5_get_F_modules.py -m *model_name* -l *sbml_level*
python P6_stripModel2Modules.py -m *model_name* -b *inf_bound* -l *sbml_level*

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, B. Teusink, and F.J. Bruggeman.

Some scripts (P5 and P6) accept additional arguments. Please read the documentation of these scripts carefully.
'

# 1. toy model
python P1_model_setup.py toy_model
python P2_find_duplicate_reactions.py -m toy_model -b 1000 -l 3
python P3_get_F_modules.py -m toy_model -l 3
python P4_SBML2H-format.py -m toy_model -b 1000 -l 3
sh ./P5_do_rational_FBA.sh toy_model
python P6_stripModel2Modules.py -m toy_model -b 1000 -l 3

# 2. iTM686.light
python P1_model_setup.py iTM686.light
python P2_find_duplicate_reactions.py -m iTM686.light -b 999999 -l 2
python P3_get_F_modules.py -m iTM686.light -l 2
python P4_SBML2H-format.py -m iTM686.light -b 999999 -l 2
sh ./P5_do_rational_FBA.sh iTM686.light
python P6_stripModel2Modules.py -m iTM686.light -b 999999 -l 2

# 3. Ecoli_iAF1260_noox.glc
python P1_model_setup.py Ecoli_iAF1260_noox.glc
python P2_find_duplicate_reactions.py -m Ecoli_iAF1260_noox.glc -b 999999 -l 2
python P3_get_F_modules.py -m Ecoli_iAF1260_noox.glc -l 2
python P4_SBML2H-format.py -m Ecoli_iAF1260_noox.glc -b 999999 -l 2
sh ./P5_do_rational_FBA.sh Ecoli_iAF1260_noox.glc
python P6_stripModel2Modules.py -m Ecoli_iAF1260_noox.glc -b 999999 -l 2
