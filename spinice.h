#ifndef SPINHEAD_H_
#define SPINHEAD_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include "randomc.h"
#include "mkl.h"

//const double schedule[] = {40, 

//double calcacorr(int N, bool* initstate, bool* spinstate);

bool flipsingle(int N, CRandomMersenne& RanGen_mersenne, int spinind, bool* spinstate, double* intmat, double T, double* energy, double* inter, double* spinstated);

//void getacorr(double* acorr, int N, bool* spinstate, double* intmat, double T, int S, int k);

double getenergy(int N, double* intmat, bool* spinstate, double* inter, double* spinstated);

double mcstep(int N, CRandomMersenne& RanGen_mersenne, double* intmat, double T, bool* spinstate, double* energy);

//void dipassign(double* dip, int m, bool ori);

void evolve(bool* state, int N, double* intmat, double temp, CRandomMersenne& RanGen_mersene);

void evolvesave(bool* timestate, int N, double* intmat, double temp, CRandomMersenne& RanGen_mersenne, int S);

#endif /* SPINHEAD_H_ */
