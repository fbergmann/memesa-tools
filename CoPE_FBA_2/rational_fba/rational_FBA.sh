#!/bin/sh
: '
Perform rational FBA on the INE model description file

Use model name as argument:  
./rational_FBA.sh model_name

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: t.r.maarleveld@cwi.nl
Last Change: July 22, 2014
'

STARTTIME=$(date +%s)

# I assume that this jar file is in the current directory, if you use
# a newer version of JLinAlg then this will need to be updated

echo $1
filepath=../data/$1/models/h-format/$1.noinf_r.ine # INE file directory
echo $filepath

JLINALGJAR=jlinalg_0.5.jar

echo $CLASSPATH | grep "$JLINALGJAR"

if [ "${?}" -ne "0" ]
	then
	echo "Couldn't find $JLINALGJAR in your Java CLASSPATH, am adding it"
	CLASSPATH=$CLASSPATH:./$JLINALGJAR
	echo "CLASSPATH is now $CLASSPATH"
	export CLASSPATH
	fi

OPTVAR=`./getObjective.sh $filepath`
echo "Variable to optimize is x$OPTVAR"

OPSTRING="x$OPTVAR"

ORIGFILE=$filepath;

cat $filepath | gawk '{if(NF>0) print $0;}' > $filepath.noblank

cat $filepath.noblank | gawk '{if($0=="end") {print $0;exit;} else print $0;}' > $filepath.nofunc

MODFILE="$filepath.nofunc"

ROWS=`cat $MODFILE | egrep '[0-9]+ +[0-9]+ +r.*' | gawk '{print $1}'`
COLS=`cat $MODFILE | egrep '[0-9]+ +[0-9]+ +r.*' | gawk '{print $2}'`

echo "The original .ine file $filepath has $ROWS rows and $COLS columns."

echo "Producing template CPLEX file in template.cplex"
cat $MODFILE | gawk -f producetemplateFULLYSYMBOLIC.gawk > template.cplex

echo "Making $OPSTRING the objective function...moving to template2.cplx"

cat template.cplex | sed "s/obj: slack/obj: $OPSTRING/" > template2.cplex

rm unnamed.sol >/dev/null 2>/dev/null

./esolver -O -L template2.cplex >/dev/null 2>/dev/null

grep 'Value' unnamed.sol
OPTIMUM=`grep 'Value' unnamed.sol | gawk '{ print $3 }'`

echo "The optimum value was according to ESOLVER was $OPTIMUM"

echo "As a floating point this is approximately..."
rm oneRational >/dev/null 2>/dev/null
echo "$OPTIMUM" > oneRational
java ToFloat oneRational

echo "output stored at", $filepath.sol

mv unnamed.sol $filepath.sol


