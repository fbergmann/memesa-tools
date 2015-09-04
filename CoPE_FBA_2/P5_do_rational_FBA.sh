#!/bin/bash
: '
Perform rational FBA (QSopt_EX solver) 

Input:
- *model_name* 

We use H-representation of the convex polyhedron

Use model_name as (first) argument:  
./P5_do_rational_FBA.sh model_name

(C) Timo R. Maarleveld, M.T. Wortel, B.G. Olivier, F.J. Bruggeman, and B. Teusink

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: February 27, 2015
'

model_name=$1

cd rational_fba
./rational_FBA.sh $model_name
