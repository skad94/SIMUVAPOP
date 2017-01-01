#ifndef Simulations_VA.h
#define Simulations_VA.h
#include <stdlib.h>
#include<stdio.h>
#include <time.h>
#include <String.h>
#include <math.h>
#include <ctype.h>
#include<string.h>

double * Simulation (int );
void Free_Simulation(double *);
void Initialise(double *, int );
double * Proc_Poisson_Simul(double , int );
void Simul_Sample (double *, int , double , double , double );
void affiche(double *, int );
double MAX(double , double );
double Frequence(double *, int , double );

void Gain_cumule(double *, int , double , double , double , double );
double * Simul_Brownien(int ,double , double );
void Markov_chain_3(int , int , double , double , double , double , double , double );
void Markov_chain(double **, int , int , int );
void M_C_Gain_cumule(double *, int ,int , double , double, double , double);
void M_C_Gain_cumule_Max(double *, int ,int , int , double , double , double , double );
double * Temps_Atteinte(int ,int , double , double , double, double);
double * Proc_O_U_Simul(int , double , double , double , double);
double * Proc_CIR_Simul(int , double , double , double , double);
double * Proc_pseudo_CEV(int , double , double , double , double);

void Fcts_Barriere_up_and_down(double *, int , double , double );
double sd_hat(double * ,int );
double * Trier(double *, int );
int TVI(double *,int ,double);

double * BS_Simul(int , double , double , double );
double * BS_Simul_milstein(int , double , double , double );
double f(double,double = 10^(-7));
double h(double);
double * Lyvat_Simul_X(int , double , double , double );
double * Lyvat_Simul_Z(int , double , double , double );
double * Lyvat_Simul_Q(int , double , double , double );
void Call_BS_MC(int , int , double , double , double, double);





#endif
