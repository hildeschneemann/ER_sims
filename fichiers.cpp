// Functions to open input and output files,
// read parameter values from input file and
// write them in output file.

#include "mutation.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;

// opens input file:

void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r"); // ERIC: Opens fichierLecture file in read mode
}


// opens output file:

void ouvrirFichierS()
{
	fichierS = fopen(fichierEcriture,"a"); // ERIC: Opens fichierEcriture file in append mode
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0
bool lireFichier(int &Nr, int &nbSr, int &NbGenr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &U_t, double &Rgr, double &Rcr, int &Rep, int &outputr, double &sdopt, double &steep)
{
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
	{
		term = true;
	}
	else
	{
		fscanf(fichierE,"%d ",&Nr);  // ERIC: "fscanf(FILE_ptr, format, address)" reads data from stream and stores at address in given format
		fscanf(fichierE,"%d ",&nbSr);
		fscanf(fichierE,"%d ",&NbGenr);
		fscanf(fichierE,"%d ",&pasr);
        fscanf(fichierE,"%lf ",&sr);
        fscanf(fichierE,"%lf ",&s_maxr);
        fscanf(fichierE,"%lf ",&Ir);
        fscanf(fichierE,"%lf ",&U_gr);
        fscanf(fichierE,"%lf ",&U_cr);
	    fscanf(fichierE,"%lf ",&U_t);
	    fscanf(fichierE,"%lf ",&Rgr);
        fscanf(fichierE,"%lf ",&Rcr);
        fscanf(fichierE,"%d ",&Rep);
        fscanf(fichierE,"%d ",&outputr);
	fscanf(fichierE,"%lf ",&sdopt);
	fscanf(fichierE,"%lf ",&steep);

		term = false;
	}
	return term;
}


// writes parameter values in output file:

void ecrireParametres(int Nv, int nbSv, int NbGenv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_t, double Rgv, double Rcv, int outputv, double sdopt, double steep)
{
	fprintf(fichierS,"\n_________________________________________\n"); // ERIC: "fprintf(FILE_ptr, string_format, var)" writes (string_format, var) to file stream
	fprintf(fichierS,"\nN = %d", Nv);
	fprintf(fichierS,", nbS = %d", nbSv);
	fprintf(fichierS,"\ngenerations = %d", NbGenv);
	fprintf(fichierS,", pas = %d", pasv);
    fprintf(fichierS,", s = %g", s);
    fprintf(fichierS,", s_max = %g", s_max);
    fprintf(fichierS,", s = %g", I);
    fprintf(fichierS,", U_g = %g", U_g);
    fprintf(fichierS,", U_c = %g", U_c);
    fprintf(fichierS,", U_t = %g", U_t);
    fprintf(fichierS,", Rg = %g", Rgv);
    fprintf(fichierS,", Rc = %g", Rcv);
    fprintf(fichierS,", output = %d", outputv);
    fprintf(fichierS,", sdopt = %g", sdopt);
    fprintf(fichierS,", steep = %g", steep);
}
