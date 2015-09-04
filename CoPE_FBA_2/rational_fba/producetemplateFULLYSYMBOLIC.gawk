

BEGIN { started = 0;
}

(index($0, "representation") > 0) && (started==0) {
started = 1;
getline;	
# we've eaten "begin"
getline;
ROWS = $1;
COLS = $2;

bigM = 1000000;

numEqu = 0;

atRow = 1;
print("\\ Will be processing " ROWS " rows and " COLS " columns");
next;
}


(started==1) {
if( NF != COLS )
	{
	print(" Not correct number of columns in row " NR "...I counted " NF);
	exit(0);
	}

for( x=1; x<=COLS; x++ )
	{
	ineq[ atRow "." x ] = "" $x;
	}
atRow++;

if( atRow == (ROWS+1) )
	{
	getline;
	if( $0 == "end" )
		{
		print("Maximize");

		printf("obj: slack");

		printf("\n");

		print("Subject To");
	
		for(x=1; x<=ROWS; x++ )
			{
			# if( x == constraint ) printf("\\ ");

			seen = 0;

			# printf("%d: ",x);

			for( c=2; c<=COLS; c++ )
				{
                                myval = ineq[ x "." c ];
                                if( c != 2 )
                                    {
				    frontsym = substr(myval,1,1);
				    if( (myval != "0") && (myval != "0.0") && (myval != "-0") && (myval != "-0.0") && (seen!=0) && (frontsym != "-")) printf(" + ");
                                    }
				if( (myval != "0") && (myval != "0.0") && (myval != "-0") && (myval != "-0.0")  )
					{
					printf(" %s x%d", myval , c-1);
					seen = 1;
					}
                                }

			if( seen == 1 )
				{
				printf(" >= ");

				fchar = substr(ineq[x "." "1"],1,1);
				
				if( fchar == "-" ) print(substr(ineq[x "." "1"], 2));
				else
				print("-" ineq[ x "." "1" ]);

	                        # print((-1)*(ineq[ x "." "1" ]));
				}
			else
				{
				print(" 0x1 >= 0");
				}

			}

		print("slack <= 1000000");
	
		print("Bounds");

		printf(" -inf <= slack <= 1000000\n");

		for(s=1; s<=(COLS-1); s++ )
			{
			printf(" -inf <= x%d <= +inf\n",s);
			}

		# print("General");
		print("End");

		exit(0);
		}
	else    {
		print("Expected an 'end' at line " atRow " , but got: " $0);
		}
	}
	else next;

}

