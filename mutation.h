// Header file: definitions of global variables, function prototypes

#ifndef MUTATION_H
#define MUTATION_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// Global variables:

#define fichierLecture "parametres.txt"     // names of input
#define fichierEcriture "resultats.txt"		// and output files

// "chr": represents a chromosome

struct ind
{
    double * gene; // allelic values at loci corresponding to fitness
    double * cis; // cistronic values at loci controlling regulation [doubles]
    double * trans; // trans values for each trans locus [CEL version!!)
    //   	Q_male =  q_x + q_y
    //      Q_female = [q1 + q2] * (n1+n2)/2

};

// Function prototypes:

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(int Nv, int nbSv, int NbGenv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_t, double Rgv, double Rgc, int outputv, double sdopt, double steep);
bool lireFichier(int &Nr, int &nbSr, int &NbGenr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &U_tr, double &Rgr, double &Rcr, int &Rep, int &outputr, double &sdopt, double &steep);
void recursion(int Nv, int nbSv, int NbGenv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_t, double Rg, double Rc, int Rep, int output, double** allAverages, double sdopt, double steep);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS);
//void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex);
double getfitness(ind &parent, double h0, double s_max, double Q0, double I, int nbS, double sdopt, double steep);
void record_averages_seq(ind * pop, double * Wtot, double * popAve, int nbSv, int Nv, double expdom, double smax, double Qopt, double I, double sdopt, double steep);

#endif
