// main() function

#include "mutation.h"
#include "MersenneTwister.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

MTRand rnd;


int main()
{

	int i, k; // ERIC: declare integer variables (see above)
	
	int nbSv = 10;
	int Nv = 2;
	double Rc = 0.1;
	double Rg = 10;
	double ML = 0;

	ind * pop = new ind [Nv]; // "pop" chr pointer refers to a chr type (see "mutation.h") with Nv dimensions; "chr" has "gene", "cis", "trans", "sex" array
	ind * temp = new ind [Nv]; // temp for saving current generation during recombination
        ind * pc; // for recombination
	
	for (i = 0; i < Nv; i++) 
	{
		//GENES
		pop[i].gene = new double [2];
        	temp[i].gene = new double [2];

		// CIS regulators
		pop[i].cis = new double [2*nbSv]; // two chromosomes
		temp[i].cis = new double [2*nbSv];

		// TRANS regulators
		pop[i].trans = new double [2*nbSv]; // two chromosomes, trans reg's  corresponding to each gene
		temp[i].trans = new double [2*nbSv];
	}
	pop[0].gene[0] = 1; // initialize all genes with 1 (WT/good) allele for Chromosome #1
	pop[0].gene[1] = 2; // same for Chromosome #2
	pop[1].gene[0] = 3;
	pop[1].gene[1] = 4;	
						
		for (k = 0; k < nbSv; k++)
		{
		    pop[0].cis[k] = 0; // initialize all cis regulators
		    pop[0].cis[nbSv+k] = 1; // same for Chromosome #2
		    pop[0].trans[k] = 0; // initialize all trans  regulators
		    pop[0].trans[nbSv+k] = 1; // same for Chromosome #2
															 }
		for (k = 0; k < nbSv; k++)
                {
                    pop[1].cis[k] = 2; // initialize all cis regulators
                    pop[1].cis[nbSv+k] = 3; // same for Chromosome #2
                    pop[1].trans[k] = 2; // initialize all trans  regulators
                    pop[1].trans[nbSv+k] = 3; // same for Chromosome #2

                 }
		
	rec(temp[0], pop[0], pop[1], Rg, Rc, ML, nbSv);

	return 0; // ERIC: C++ Convention (void main() prohibited); indicates program successful
}
