/*
 * Assemble Hamiltonian matrix (tight-binding + spin-orbit interaction)
 */
#include "slater-koster.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdio.h>


#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#define SPIN 2

double cos_dir(double, double, double, double );

gsl_matrix_complex * hamiltonian(double **, int N);

#endif // hamiltonian
