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

//const double schedule[] = {40};

const int cellsize = 8;

const int N = 16*cellsize*cellsize*cellsize;

bool flipsingle(CRandomMersenne& RanGen_mersenne, int spinind, bool* spinstate, double* intmat, double T, double* energy, double* inter, double* spinstated);

double getenergy(double* intmat, bool* spinstate, double* inter, double* spinstated);

double mcstep(CRandomMersenne& RanGen_mersenne, double* intmat, double T, bool* spinstate, double* energy, int k);

void evolve(bool* state, double* intmat, double temp, CRandomMersenne& RanGen_mersene,int label);

void evolvesave(bool* timestate, double* intmat, double temp, CRandomMersenne& RanGen_mersenne, int S, int label);

#endif /* SPINHEAD_H_ */
