#include "mutation.h"
#include <cmath>
using namespace std;

#define PIE 3.14159265359
/*void record_output(ind * pop, double * Wtot, double ** measures, double * popAve, int nbSv, int Nmales, int Nv, double expdom)
{
    // measures from population: for each locus: sbarX, sbarY, hY, attractX, cisx,cisy,trans

    double WbarMales, WbarFemales, y, transCel, c1, c2;
    int i, j;
    int Nfemales = Nv - Nmales;
    int twoNfemale = 2*(Nv - Nmales);

    for (j = 0; j < nbSv; j++)
        for (i = 0; i < 7; i++)
            measures[j][i] = 0;

    WbarMales = 0;
    for (i = 0; i < Nmales; i++)
    {
        WbarMales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transCel = (((pop[i].trans[j] > 0) ? pop[i].trans[j] : 0)+((pop[i].trans[nbSv+j] > 0) ? pop[i].trans[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            measures[j][0] += 1 - pop[i].gene[nbSv+j];
            measures[j][1] += 1 - pop[i].gene[j];
            y = c1 / (c1 + c2);
            measures[j][2] += pow(y, expdom);
            measures[j][4] += c2;
            measures[j][5] += c1;
            measures[j][6] += transCel;

        }
    }
    WbarMales /= Nmales;

    WbarFemales = 0;
    for (i = Nmales; i < Nv; i++)
    {
        WbarFemales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transCel = (((pop[i].trans[j] > 0) ? pop[i].trans[j] : 0)+((pop[i].trans[nbSv+j] > 0) ? pop[i].trans[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            measures[j][0] += 2 - pop[i].gene[j] - pop[i].gene[nbSv+j];
            measures[j][3] += (c1 + c2) * transCel;
            measures[j][4] += (c1 + c2) ;
            measures[j][6] +=  transCel ;

        }

    }
    for (j = 0; j < nbSv; j++)
    {
        measures[j][0] /= (Nmales + twoNfemale); // sbarX
        measures[j][1] /= Nmales; // sbarY
        measures[j][2] /= Nmales; // hY
        measures[j][3] /= twoNfemale; // attractX
        measures[j][4] /= (Nmales + twoNfemale); // cisX
        measures[j][5] /= Nmales; // cisY
        measures[j][6] /= (Nmales+Nfemales); // transcel

    }
    WbarFemales /= Nfemales;

    popAve[0] = WbarMales;
    popAve[1] = WbarFemales;

}

void record_averages(ind * pop, double * Wtot, double * popAve, int nbSv, int Nv, double expdom, double smax)
{
    // measures from population: Wbarmale, WbarFemale, sbarXbar, sbarYbar, hYbar, attractXbar, cisxbar,cisybar,transbarbar,
    // number of loci half-dead, dead, half-silenced, silenced on Y

    double Wbar, Celbar, y, transCel, c1, c2;
    int i, j;
    double halfminW = (1-smax/2.0);
    double minW = 1-smax;
    int cmptyM = 0;


    for (i = 0; i < 13; i++)
        popAve[i] = 0;

    Wbar = 0;
    Celbar = 0;
    for (i = 0; i < Nv; i++)
    {
        Wbar += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transCel = (((pop[i].trans[j] > 0) ? pop[i].trans[j] : 0)+((pop[i].trans[nbSv+j] > 0) ? pop[i].trans[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            popAve[3] += 1 - pop[i].gene[j];       //sbarybar
			
            if (c1 + c2 > 0)
            {
                y = c1 / (c1 + c2);
                popAve[4] += pow(y, expdom);    //hybar
                cmptyM++;
            }
            popAve[6] += c2;  //cisxbar
            if (pop[i].gene[j] < halfminW)
                popAve[9] += 1;
            if (pop[i].gene[j] == minW)
                popAve[10] += 1;
            if (y < 0.25)
                popAve[11] += 1;
            if (y < 0.01)
                popAve[12] += 1;
        }
    }
    Wbar /= Nv;


    popAve[0] = Wbar;

    popAve[2] /= (nbSv *(Nv + twoNfemale)); // sbarXbar;
    popAve[3] /= (nbSv * Nmales); // sbarYbar ;
    popAve[4] /= cmptyM;  //hY
    popAve[5] /= (nbSv * twoNfemale); //attractXbar
    popAve[6] /= (nbSv *(Nmales + twoNfemale)); //cisxbar
    popAve[7] /= (nbSv *Nmales); // cisybar
    popAve[8] /= (nbSv * Nv); // transbar
    popAve[9] /= (nbSv * Nmales);
    popAve[10] /= (nbSv * Nmales);
    popAve[11] /= (nbSv * Nmales);
    popAve[12] /= (nbSv * Nmales);

}

*/
void record_averages_seq(ind * pop, double * Wtot, double * popAve, int nbSv, int Nv, double expdom, double smax, double Qopt, double I, double sdopt, double steep)
{
	double Wbar, e1, e2, eopt;
	int i,j, hom, het, mm11, mm12, mm21, mm22, allele1, allele2, matches;
	hom = 0;
	het = 0;
	Wbar = 0.0;
	matches = 0;
	
	for (i = 0; i < Nv; i++)
	{
        	Wbar += Wtot[i];
		if (pop[i].gene[0] == pop[i].gene[1])
			hom += 1;
		else
			het += 1;
			
		mm11 = 0;
		mm12 = 0;
		mm21 = 0;
		mm22 = 0;

        	for (j = 0; j < nbSv; j++)
        	{
			allele1 = j;
			allele2 = nbSv + j;
			//record genotypes individuals
			//fout2 << pop[i].cis[allele1] << pop[i].cis[allele2] << pop[i].trans[allele1] << pop[i].trans[allele2];	

			mm11 += (pop[i].cis[allele1] == pop[i].trans[allele1] ? 0 : 1);
			mm12 += (pop[i].cis[allele1] == pop[i].trans[allele2] ? 0 : 1);
			mm21 += (pop[i].cis[allele2] == pop[i].trans[allele1] ? 0 : 1);
			mm22 += (pop[i].cis[allele2] == pop[i].trans[allele2] ? 0 : 1);
		}
		eopt = Qopt + sqrt(1.0/(2.0*I)) * sdopt / (1/(steep * sqrt(2 * PIE)));
		e1 =  	(1/(steep * sqrt(2 * PIE))) * (exp(-mm11*mm11/(2 * (steep*steep)))) + 
			(1/(steep * sqrt(2 * PIE))) * (exp(-mm12*mm12/(2 * (steep*steep))));
		e2 =  	(1/(steep * sqrt(2 * PIE))) * (exp(-mm21*mm21/(2 * (steep*steep)))) + 
			(1/(steep * sqrt(2 * PIE))) * (exp(-mm22*mm22/(2 * (steep*steep))));
		e1 *= eopt;
		e2 *= eopt;	
		matches += (mm11 + mm12 + mm21 + mm22);
 
	}
	//fout2 << endl;
	Wbar /= Nv;
	popAve[0] = Wbar; // mean fitness
	popAve[1] = double(matches) / double(Nv*4.0); // mean match
	popAve[2] = double(hom) / double(Nv); // frequency homozygotes at gene locus
	popAve[3] = double(het) / double(Nv); // frequency heterozygotes at gene locus
	
}
