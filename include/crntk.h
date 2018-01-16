/*******************************************************************************
 * CRNTK: The Chemical Reaction Network Toolkit                                *
 * Copyright (C) 2018 Jason Dark, email@jkdark.com                             *
 *                                                                             *
 * Redistribution and use in source and binary forms, with or without          *
 * modification, are permitted provided that the following conditions are met: *
 *                                                                             *
 * 1. Redistributions of source code must retain the above copyright notice,   *
 *    this list of conditions and the following disclaimer.                    *
 *                                                                             *
 * 2. Redistributions in binary form must reproduce the above copyright notice,*
 *    this list of conditions and the following disclaimer in the documentation*
 *    and/or other materials provided with the distribution.                   *
 *                                                                             *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" *
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   *
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         *
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        *
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    *
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  *
 * POSSIBILITY OF SUCH DAMAGE.                                                 *
 ******************************************************************************/


#ifndef CRNTK_H
#define CRNTK_H

#include <stdbool.h>
#include <stddef.h> // size_t
#include <stdlib.h> // calloc
#include <string.h> // memset, memcpy


// CRNTK provides several common kinetic schemes, and users may implement their provide their own (e.g. a hill function)
// For the user-implemented case, the 3rd argument is passed a runtime-fixed array of real numbers that can be used however
double crntk_kinetics_mass_action(size_t, size_t, const double*);
double crntk_kinetics_heaviside(size_t, size_t, const double*);
double crntk_kinetics_fsp(size_t, size_t, const double*);

// The crntk_reactant represents a distinct reactant with specified kinetics
typedef struct  {
    double (*affinity)(size_t,size_t,const double*);
    double *data;
} crntk_reactant;

// The crntk_reaction associates a rate with two complexes
typedef struct {
    double rate;
    size_t* lhs;
    size_t* rhs;
} crntk_reaction;

// The crntk_crn wraps every piece of data together, then serves as the argument to pass to the operations
typedef struct {
    // the "dimensions" of the CRN
    size_t n_reactants;
    size_t n_complexes;
    size_t n_reactions;
    size_t n_states;
    size_t n_constraints;

    crntk_reactant *reactants;
    size_t         *complexes;
    crntk_reaction *reactions;
    size_t         *states;
    size_t         *constraints;
    size_t         *constraint_values;

    // a representation of the lattice
	size_t **table;
} crntk_crn;

// used as an internal method, can be used externally. Given an array [n1,n2,n3,...], gives the integer
// offset corresponding to that state. Not defined if the state does not satisfy the conservation laws
size_t crntk_index_of(const crntk_crn *, const size_t *);


// Finally, these 2 are the matrix methods intended for normal user consumption.
// Suppose the CME reads p'=Ap:
void crntk_id_apply(        const crntk_crn* crn, const double *x, double *y); // Computes y=Ax
void crntk_tr_apply(        const crntk_crn* crn, const double *x, double *y); // Computes y=xA
// Let tr(A)=L+D+U. To compute absorption times, we might need to precondition. Here are 4 preconditioners:
void crntk_tr_sor_forward(  const crntk_crn* crn, double omega, double *x);    // Computes x=inv(w*L-D)x
void crntk_tr_sor_backward( const crntk_crn* crn, double omega, double *x);    // Computes x=inv(w*U-D)x
void crntk_tr_ssor_forward( const crntk_crn* crn, double omega, double *x);    // Computes x=inv(w*U-D)inv(D)inv(w*L-D)x
void crntk_tr_ssor_backward(const crntk_crn* crn, double omega, double *x);    // Computes x=inv(w*L-D)inv(D)inv(w*U-D)x

// These are to initialize and destroy allocated memory. You should always call these.
bool crntk_init(crntk_crn *lattice);
void crntk_destroy(crntk_crn *lattice);

#endif
