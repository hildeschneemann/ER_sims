#include "mutation.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

// rec function: generates recombinant genome segment "res"
// from parental genome segments "c1" and "c2"
// nbCo is the number of cross-overs in the genome segment,
// nS the number of selected loci, nU nb of loci affecting U


/*void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex)
{
    vector<double> posCoG;
    int i, j, locus, nbCo, ChrInit, strand;

    // paternally inherited gamete
    double Rsex = Rg; // NEW linkage between the first cis reg and the sex locus. Only matters in male where the sex locus is heterozygous. Set to Rg to avoid adding a parameter
    nbCo = poisdev(ML+Rsex); // NEW number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    if (off_sex == 0) // if offspring is a male, it inherits the sex-determining locus from the proto-Y
        ChrInit = 0;
    else
        ChrInit = 1;

    // CIS REG:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while (locus * Rg + Rsex < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        {
            offspring.cis[locus] = father.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    {
        offspring.cis[locus] = father.cis[strand + locus];
        locus++;
    }

    // GENES:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus * Rg + Rc + Rsex) < posCoG[j]) // NEW note that position of gene i is i*Rg + Rc
        {
            offspring.gene[locus] = father.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[locus] = father.gene[strand + locus];
        locus++;
    }

    // maternally inherited gamete

    nbCo = poisdev(ML); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    ChrInit = rnd.randInt(1); // which X chromosome contributes at position 0

    // CIS REG:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while (locus * Rg < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        {
            offspring.cis[nS+locus] = mother.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    {
        offspring.cis[nS+locus] = mother.cis[strand + locus];
        locus++;
    }

    // GENES:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus * Rg + Rc) < posCoG[j]) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[nS+locus] = mother.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[nS+locus] = mother.gene[strand + locus];
        locus++;
    }

    // TRANS REGULATORS

    for (i = 0; i < nS; i++)
    {
        // paternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i] = father.trans[ChrInit*nS+i];

        // maternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i+nS] = mother.trans[ChrInit*nS+i];
    }
}

*/

void rec(ind &offspring, ind &parent1, ind &parent2, double Rg, double Rc, double ML, int nS)
{
    	vector<double> posCoG;
    	int i, j, locus, nbCo, ChrInit, strand;
	double pos;

/*
	//print genotypes of parents
	cout << "genotype parent 1: \n";
	for (i = 0; i < nS; i++)
	{
		cout << parent1.cis[i];
	}
	cout << " " << parent1.gene[0] << endl;

	for (i = 0; i < nS; i++)
        {
                cout << parent1.cis[nS + i];
        }
        cout << " " << parent1.gene[1] << endl;

	cout << "genotype parent 2: \n";
        for (i = 0; i < nS; i++)
        {
                cout << parent2.cis[i];
        }
        cout << " " << parent2.gene[0] << endl;

        for (i = 0; i < nS; i++)
        {
                cout << parent2.cis[nS + i];
        }
        cout << " " << parent2.gene[1] << endl;

	cout << endl;
*/
	//PARENT 1 recombination
	
	nbCo = poisdev(Rc*nS); // number of cross-overs
        //cout << "PARENT 1 nbCo: " << nbCo << ", positions:";

    	for (j = 0; j < nbCo; j++)
        {
                pos = rnd.randExc(nS);
                //cout << " " << pos;
                posCoG.push_back(pos);
        }
        //cout << endl;
    	sort(posCoG.begin(), posCoG.end());

    	ChrInit = rnd.randInt(1);
	//cout << "ChrInit: " << ChrInit;
	//cout << " strand: ";
    	locus = 0; // ERIC: initialize locus again
    	for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    	{
        	strand = ((j+ChrInit)%2)*nS;
        	while (locus < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        	{
			//cout <<  strand;
            		offspring.cis[locus] = parent1.cis[strand + locus];
            		locus++;
        	}
    	}
    	strand = ((nbCo+ChrInit)%2)*nS;
    	while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    	{
        	offspring.cis[locus] = parent1.cis[strand + locus];
        	locus++;
		//cout << strand;
    	}
	//cout << endl;
    	nbCo = poisdev((Rc*Rg));
    	offspring.gene[0] = parent1.gene[((nbCo+strand/nS)%2)];
	
	//cout << "PARENT1 nbCo gene: " << nbCo << endl;
	//cout << endl;
	posCoG.clear();
    	// PARENT2 recombination between chromosomes parent 2
    	nbCo = poisdev((Rc*nS)); // number of cross-overs
    	for (j = 0; j < nbCo; j++)
        	posCoG.push_back(rnd.randExc(nS)); // position of cross-over j
    	sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

	/*
	cout << "PARENT2 nbCo: " << nbCo << ", positions: ";
        for (j = 0; j < nbCo; j++)
        {
                cout << " " << posCoG[j];
        }
        cout << endl;
	*/
    	ChrInit = rnd.randInt(1);

	//cout << "ChrInit: " << ChrInit << " strand: ";
    	locus = 0; // ERIC: initialize locus again
    	for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    	{
        	strand = ((j+ChrInit)%2)*nS;
        	while (locus  < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        	{
			//cout << strand;
            		offspring.cis[nS+locus] = parent2.cis[strand + locus];
            		locus++;
        	}
    	}
    	strand = ((nbCo+ChrInit)%2)*nS;
    	while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    	{
		//cout << strand;
        	offspring.cis[nS+locus] = parent2.cis[strand + locus];
        	locus++;
    	}
	//cout << endl;
    	nbCo = poisdev((Rc*Rg));
	
 	//cout << "PARENT2 nbCo gene: " << nbCo << endl;
	//cout << endl;
    	offspring.gene[1] = parent2.gene[((nbCo+strand/nS)%2)];

	/*
	cout << "genotype offspring: \n";
        for (i = 0; i < nS; i++)
        {
                cout << offspring.cis[i];
        }
        cout << " " << offspring.gene[0] << endl;

        for (i = 0; i < nS; i++)
        {
                cout << offspring.cis[nS + i];
        }
        cout << " " << offspring.gene[1] << endl;


	//cout << "finished CIS & GENE recombination\n";
	cout << endl;
    	*/

    	// TRANS REGULATORS
    	// recombination between chromosomes parent 1 TRANS region
	posCoG.clear();
 	nbCo = poisdev((Rc*nS)); // number of cross-overs
    	for (j = 0; j < nbCo; j++)
        	posCoG.push_back(rnd.randExc(nS)); // position of cross-over j
    	sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last
	
	/*
	cout << "PARENT 1 TRANS nbCo: " << nbCo << ", positions: ";
	for (j = 0; j < nbCo; j++)
	{
		cout << " " << posCoG[j];
	}
	cout << endl;
	*/
    	ChrInit = rnd.randInt(1);
	//cout << "ChrInit: " << ChrInit << endl;

    	locus = 0; // ERIC: initialize locus again
    	for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    	{
        	strand = ((j+ChrInit)%2)*nS;
        	while (locus  < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        	{
			//cout << strand;
            		offspring.trans[locus] = parent1.trans[strand + locus];
            		locus++;
        	}
    	}	
    	strand = ((nbCo+ChrInit)%2)*nS;
    	while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    	{
		//cout << strand;
        	offspring.trans[locus] = parent1.trans[strand + locus];
        	locus++;
    	}
	//cout << endl;

	posCoG.clear();

    	//recombination between chromosomes parent 2 TRANS regulator
    	posCoG.clear();
	nbCo = poisdev((Rc*nS)); // number of cross-overs
    	for (j = 0; j < nbCo; j++)
        	posCoG.push_back(rnd.randExc(nS)); // position of cross-over j
    	sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    	ChrInit = rnd.randInt(1);
	/*cout << "PARENT 2 TRANS nbCo: " << nbCo << ", positions: ";
        for (j = 0; j < nbCo; j++)
        {
                cout << " " << posCoG[j];
        }
        cout << endl;
	*/
    	ChrInit = rnd.randInt(1);
        //cout << "ChrInit: " << ChrInit << endl;


    	locus = 0; // ERIC: initialize locus again
    	for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    	{
        	strand = ((j+ChrInit)%2)*nS;
        	while (locus  < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        	{
			//cout << strand;
            		offspring.trans[nS+locus] = parent2.trans[strand + locus];
            		locus++;
        	}
    	}
    	strand = ((nbCo+ChrInit)%2)*nS;
    	while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    	{
		//cout << strand;
        	offspring.trans[nS+locus] = parent2.trans[strand + locus];
        	locus++;
    	}
    	//cout << endl;

/*
	// RECORD TRANS GENOTYPES PARENTS AND OFFSPRING
	cout << endl;
	cout << "genotype parent 1 trans: \n";
        for (i = 0; i < nS; i++)
        {
                cout << parent1.trans[i];
        }
	cout << endl;

        for (i = 0; i < nS; i++)
        {
                cout << parent1.trans[nS + i];
        }

	cout << endl;
        cout << "genotype parent 2 trans: \n";
        for (i = 0; i < nS; i++)
        {
                cout << parent2.trans[i];
        }
	
	cout << endl;
        for (i = 0; i < nS; i++)
        {
                cout << parent2.trans[nS + i];
        }

        cout << endl;
	cout << endl;
	//cout << "finished trans recombination\n";
	cout << "genotype offspring trans: \n";
        for (i = 0; i < nS; i++)
        {
                cout << offspring.trans[i];
        }
	cout << endl;

        for (i = 0; i < nS; i++)
        {
                cout << offspring.trans[nS + i];
        }
	cout << endl;
*/
}
