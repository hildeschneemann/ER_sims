#include "mutation.h"
#include <cmath>
using namespace std;

#define PIE 3.14159265359
extern MTRand rnd;
/*double Wmale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j, allele1, allele2;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus,explvscaled;


    for (j = 0; j < nbS; j++)
    {
        allele1 = j;
        allele2 = nbS+j;

        e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);  // cis_y
        e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0) ;  // (m1+m2)/2 * cis_x

        explv = e1 + e2; //  cis_x + cis_y
        //stabilizing selection fitness effect
        //                   w = w * exp(-I*pow((explv-Qopt),2)); // OLD

        if (explv == 0)
            w *= 1 - smax;

        else
        {
            //effectLocus = 1-smax*(1-exp(-I*pow((explv-Qopt),2))); // OLD Gaussian version
            //effectLocus = 1-smax*(1-explv*exp(1-explv/Qopt)/Qopt);      // OLD gamma version
            if (explv < Qopt)
                explvscaled= Qopt/explv - 1;
            else
                explvscaled= explv/Qopt -1 ;   //This trick scales the fitness effect with the ratio of expression compared to Qopt, in one direction or another, and insures that fitness is always 1-smax with 0 expression

            effectLocus = 1-smax*(1-exp(-I*explvscaled*explvscaled));      // NEW



            if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e1/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[nbS+j] + h * (parent.gene[j] - parent.gene[nbS+j]); <-- additive fitness
                effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2]));
            }
            else // if first allele is fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e2/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[j] + h * (parent.gene[nbS+j] - parent.gene[j]); <-- additive fitness
                effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1]));
            }
            if (effectLocus < 1 - smax)
                w *= 1 - smax;
            else
                w *= effectLocus;
        }
    }
    return w;
}


double Wfemale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus,explvscaled,transCel;
    int allele1, allele2;


   // double transCel = (((parent.trans[3] > 0) ? parent.trans[3] : 0) + ((parent.trans[8] > 0) ? parent.trans[8] : 0))/2;



        for (j = 0; j < nbS; j++)
        {
            allele1 = j;
            allele2 = nbS+j; // to avoid doing the sum multiple times
            transCel = (((parent.trans[allele1] > 0) ? parent.trans[allele1] : 0) + ((parent.trans[allele2] > 0) ? parent.trans[allele2] : 0))/2;

            // (n1+n2)/2 * [(o1+o2)/2*cis_x1 + cis_x2]
            e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);
            e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0);
            explv = (e1 + e2) * transCel ;

            if (explv == 0)
                w *= 1 - smax;

            else
            {
                //stabilizing selection fitness effect
                //                    w = w * exp(-I*pow((explv-Qopt),2)); // OLD
                //effectLocus = 1-smax*(1-exp(-I*pow((explv-Qopt),2))); // OLD
                //    effectLocus = 1-smax*(1-explv*exp(1-explv/Qopt)/Qopt);      // OLD

                if (explv < Qopt)
                    explvscaled= Qopt/explv - 1;
                else
                    explvscaled= explv/Qopt -1 ;   //This trick scales the fitness effect with the ratio of expression compared to Qopt, in one direction or another, and insures that fitness is always 1-smax with 0 expression

                effectLocus = 1-smax*(1-exp(-I*explvscaled*explvscaled));      // NEW


                if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e1/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2])); // product fitnesses
                }
                else // if first allele fittest
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e2/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1])); // product fitnesses
                }
                if (effectLocus < 1 - smax)
                    w *= 1 - smax;
                else
                    w *= effectLocus;
            }
        }



    return w;
}
*/
double getfitness(ind &parent, double expdom, double smax, double Qopt, double I, int nbS, double sdopt, double steep)
{
	
	int j, allele1, allele2, mm11, mm12, mm21, mm22;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus,explvscaled, eopt;
    mm11 = 0;
    mm12 = 0;
    mm21 = 0;
    mm22 = 0;

	for (j = 0; j < nbS; j++)
        {
			allele1 = j;
			allele2 = nbS + j;
			
			//compute match between cis allele 1 and both transcription factors
			mm11 += (parent.cis[allele1] == parent.trans[allele1] ? 0 : 1);
			mm12 += (parent.cis[allele1] == parent.trans[allele2] ? 0 : 1);
			mm21 += (parent.cis[allele2] == parent.trans[allele1] ? 0 : 1);
			mm22 += (parent.cis[allele2] == parent.trans[allele2] ? 0 : 1);
		}

	eopt = Qopt + sdopt * sqrt(1/(2 *I));
	e1 =  eopt * exp(- steep * mm11*mm11) + eopt * exp(- steep * mm12*mm12);
	e2 =  eopt * exp(- steep * mm21*mm21) + eopt * exp(- steep * mm22*mm22);
	explv = e1 + e2;
	
	if (explv == 0)
		w = 1 - smax;

	else
	{
		if (explv < Qopt)
			explvscaled= Qopt/explv - 1;
		else
			explvscaled= explv/Qopt -1 ;   //This trick scales the fitness effect with the ratio of expression compared to Qopt, in one direction or another, and insures that fitness is always 1-smax with 0 expression

		effectLocus = 1-smax*(1-exp(-I*explvscaled*explvscaled));      // NEW



		if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
		{
			// h is the dominance coefficient, H0 = 0.25
			h = pow((e1/(e1+e2)), expdom);

			// fitness = w1 + h*(w2-w1)
			effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2]));
		}
		else // if first allele is fittest
		{
			// h is the dominance coefficient, H0 = 0.25
			h = pow((e2/(e1+e2)), expdom);

			// fitness = w1 + h*(w2-w1)
			effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1]));
		}
		if (effectLocus < 1 - smax)
			w = 1 - smax;
		else
			w = effectLocus;
	}
	return w;
}


