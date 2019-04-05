// main() function

#include "mutation.h"
#include "MersenneTwister.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

// ERIC: The main() function calls on "MersenneTwister.h" for random seeding and "mutation.h" for the modelling


// random number generator (MTRand class defined in MersenneTwister.h)

MTRand rnd;

// pointers to input and output files:

FILE * fichierE; // ERIC: fichierE is a FILE type pointer
FILE * fichierS; // ERIC: fichierS is a FILE type pointer

int main()
{
	// definitions of variables:

	// START ERIC:
	// Nt - population size
	// sig - standard deviation of effects on phenotype
	// om2 - selection strength
	// uInit - initial mutation rate
	// nbS - number of loci that affect phenotype per genome
	// L - genome length (number of crossovers per meiotic event)
	// NbGen - number of generations
	// pas - time steps between measurements
	// alpha - rate of change of best phenotype
	// nbU - number of mutation modifiers
	// NbGenPreli - number of generations before mutation
	// sigU - standard deviation of effects on mutation rate
	// c - cost of reproduction
	// NEW: s - mean mutational effect
	// END ERIC

	int repet,i, j, Nt, nbS, NbGen, pas, Rep, output, p; // ERIC: declare integer variables (see above)
	double  s, s_max, I, U_g, U_c, U_t, Rg, Rc, sdopt, steep; // ERIC: declare double variables (see above)

	// opens input and output files:

    bool fin;
    ouvrirFichierE();
    ouvrirFichierS();
    fin = false;
    ofstream fout;

    // reads parameter values from input file, and writes them in output file:

    do
    {
        fin = lireFichier(Nt, nbS, NbGen, pas, s, s_max, I, U_g, U_c, U_t, Rg, Rc, Rep, output, sdopt, steep); // ERIC: see "fichiers.cpp"; reads parameters from input file



        if (!fin)
        {
            double** allAverages = new double *[(NbGen / pas) + 1];
            for(i = 0; i < (NbGen / pas) + 1; i++)
                allAverages[i] = new double[13];

            for(i = 0; i < (NbGen / pas) + 1; i++)
                for (j = 0; j < 13; j++)
                    allAverages[i][j] = 0;

            ecrireParametres(Nt, nbS, NbGen, pas, s, s_max, I, U_g, U_c, U_t, Rg, Rc, output, sdopt, steep); // ERIC: see "fichiers.cpp"; writes to the output file

            // runs the simulation:
            // ERIC: "recursion" is defined in "mutations.h" and is explicit in "recursion_bis.cpp"
            // ERIC: This appears to be the meat of the simulation
            //cout << "Nt: " << Nt << "sig: " << sig << "nbS: " << nbS << NbGen << pas << alpha << s << U_g << U_c << U_t;
            for (repet = 0; repet < Rep; repet++)
            {
                recursion(Nt, nbS, NbGen, pas, s, s_max, I, U_g, U_c, U_t, Rg, Rc, repet, output, allAverages, sdopt, steep); // ERIC: see "recursion_bis.cpp"
                /*char nomFichier[256];
                stringstream nomF;
                nomF << "summary_N" << Nt << "_nbS" << nbS << "_s" << s << "_I" << I << "_Ug" << U_g << "_Uc" << U_c << "_Ut" << U_t << "_Rg" << Rg << "_Rc" << Rc <<".txt"; // results file naming convention
                nomF >> nomFichier; // Writing "nomF" to nomFichier file
                fout.open(nomFichier);
                fout<<repet+1<< endl;
                for(i = 0; i < (NbGen / pas) + 1; i++)
                {
                    fout << i*pas << " ";
		    for (j = 0; j < 13; j++)
                        fout << allAverages[i][j] << " ";
                    fout << endl;
                }
                fout.close();*/
            }

            for(i = 0; i < (NbGen/pas + 1); i++)
                delete [] allAverages[i];
            delete [] allAverages;
        }
    } while (!fin);

	// closes files:

	fclose(fichierE);
	fclose(fichierS);

	return 0; // ERIC: C++ Convention (void main() prohibited); indicates program successful
}
