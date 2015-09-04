import org.jlinalg.LinSysSolver;
import org.jlinalg.Matrix;
import org.jlinalg.Rational;
import org.jlinalg.FieldElement;
import org.jlinalg.Vector;

import java.util.*;
import java.io.*;
import java.math.BigInteger;

/**
 * Example that computes a solution of a linear equation system for the domain
 * of rational numbers.
 * 
 * @author Andreas Keilhauer, Georg Thimm
 */
public class ToFloat
{
        public static boolean DEBUG = false;

	/**
	 * start the demonstration
	 * 
	 * @param argv
	 *            is ignored
	 */
	public static void main(String[] argv)
	{
	String fileName = argv[0];


                        try     {
                                FileReader fr = new FileReader(fileName);
                                BufferedReader br = new BufferedReader(fr);

                                String record = new String();
                                while((record=br.readLine()) != null)
                                        {
                                        if( record.startsWith("//") || record.startsWith("*") || record.startsWith("begin") || 
record.startsWith("H-representation") || record.startsWith)
						{
						System.out.println(record);
						continue;
						}

                                        String[] tripData = record.split(" ");

					for( int x=0; x<tripData.length; x++ )
						{
						String numerator = null;
						String denominator = null;

						String ratString = tripData[x];

						int locbar = ratString.indexOf('/');
						if( locbar == -1 )
							{
							numerator = ratString;
							denominator = "1";
							}
						else
							{
							numerator = ratString.substring( 0, locbar );
							denominator = ratString.substring( locbar+1 );
							}

						//! System.out.println("Got rational ["+numerator+"]["+denominator+"]");

						BigInteger bnum = new BigInteger(numerator);
						BigInteger dnum = new BigInteger(denominator);

						Rational myrat = new Rational(bnum, dnum);

						if( x !=0 ) System.out.print(" ");
						System.out.print(myrat.doubleValue());
						}
					System.out.println();
					}

                                }
                        catch (IOException e)
                                {
                                // catch possible io errors from readLine()
                                System.out.println("Problem reading file "+fileName);
                                System.exit(0);
                                }

		

	}
}

