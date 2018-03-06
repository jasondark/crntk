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

#include <stdbool.h> // bool
#include <stddef.h>  // size_t

/// @brief an opaque type used to represent a chemical reaction network
typedef struct crntk_crn* crntk;

#ifdef __cplusplus
extern "C" {
#endif

/* The next group of functions defines the interface to creating a CRN */

/// @brief initialize the memory for a chemical reaction network (CRN)
/// @param[inout] crn           a pointer to the memory
/// @param[in]    n_reactants   the number of unique reactants (species) in the CRN
/// @param[in]    n_complexes   the number of distinct reaction complexes in the CRN
/// @param[in]    n_reactions   the number of reactions in the CRN
/// @param[in]    n_constraints the number of conservation laws embedded in the CRN
/// @return                     a boolean success parameter
bool crntk_init(crntk* crn, const size_t n_reactants, const size_t n_complexes, const size_t n_reactions, const size_t n_constraints);

/// @brief define the next reaction in the chemical reaction network (CRN)
/// @param[inout] crn      an opaque pointer
/// @param[in]    affinity the kinetics associated with the reactant (e.g. &crntk_kinetics_mass_action)
/// @param[in]    data     arbitrary data to be passed to @p affinity when evaluated
/// @return                a unique number representing the reactant
size_t crntk_add_reactant(crntk crn, double (*affinity)(size_t,size_t,const double*), double *data);

/// @brief define the next complex in the chemical reaction network (CRN)
/// @param[inout] crn  an opaque pointer
/// @param[in]    ...  the multipliers of each (ordered) reactant in this complex (e.g. for reactants A, B, C, the complex 2B+C would be represented by 0, 2, 1)
/// @return            a unique number representing the complex -- useful when used in conjunction with `crntk_add_reaction`
size_t crntk_add_complex(crntk crn, ...);

/// @brief define the next reaction in the chemical reaction network (CRN)
/// @param[inout] crn  an opaque pointer
/// @param[in]    rate the stochastic rate constant (units of 1/time)
/// @param[in]    lhs  the offset of the source complex
/// @param[in]    rhs  the offset of the destination complex
/// @return            a unique number representing the reaction -- useful when redefining a reaction's rate constant
size_t crntk_add_reaction(crntk crn, double rate, size_t lhs, size_t rhs);

/// @brief modify an already-defined reaction to have the specified rate
/// @param[inout] crn  an opaque pointer
/// @param[in]    rxn  the reaction id to modify
/// @param[in]    rate the new rate constant
void crntk_set_reaction_rate(crntk crn, size_t rxn, double rate);

/// @brief define the next linear constraint in the chemical reaction network (CRN)
/// @param[inout] crn   an opaque pointer
/// @param[in]    value the value that the inner product must obtain
/// @param[in]    ...   the multipliers of each (ordered) reactant in the conservation law
/// @return             a unique number representing the constraint
size_t crntk_add_constraint(crntk crn, size_t value, ...);

/// @brief finalize the chemical reaction network, after specifying each of its components
/// @param[inout] crn an opaque pointer
/// @return           a boolean success parameter
bool crntk_finalize(crntk crn);

/// @brief frees any memory that was allocated by crntk functions and owned by crn
/// @param[inout] crn an opaque pointer
void crntk_destroy(crntk crn);




/* the next 4 functions are helper functions for reading-out quantities of interest */

/// @brief computes the array offset of state @p n
/// @param[in] crn an opaque pointer
/// @param[in] n   the copy numbers of each reactant
/// @return        the array offset
size_t crntk_index_of(const crntk crn, const size_t* n);

/// @brief returns the number of discrete copy number configurations satisfying the conservation laws
/// @param[in] crn an opaque pointer
/// @return        the number of dimensions
size_t crntk_dim(const crntk crn);

/// @brief returns a pointer to the ith discrete state
/// @param[in] crn an opaque pointer
/// @param[in] i   the integer offset
/// @return        the pointer to the state
const size_t* crntk_state(const crntk crn, const size_t i);

/// @brief fills a vector with the diagonal elements of the matrix
/// @param[in]    crn an opaque pointer
/// @param[inout] d   the vector to fill
void crntk_diag(const crntk crn, double* d);




/* here are some builtin kinetic schemes to be used by the user -- note that *
 * custom kinetics are fairly straight-forward to implement                  */

/// @brief standard mass-action kinetics (e.g. 2X -> ... has rate n*(n-1)/2))
/// @param[in] x the current copy number
/// @param[in] n the multiplier
/// @param[in] data (ignored)
double crntk_kinetics_mass_action(size_t x, size_t n, const double* data);

/// @brief heaviside kinetics (e.g. 2X -> ... has rate 1 if n >=2 or 0 otherwise)
/// @param[in] x the current copy number
/// @param[in] n the multiplier
/// @param[in] data (ignored)
double crntk_kinetics_heaviside(size_t x, size_t n, const double* data);

/// @brief constant kinetics (e.g. 2X -> ... always has rate 1)
/// @param[in] x the current copy number
/// @param[in] n the multiplier
/// @param[in] data (ignored)
double crntk_kinetics_fsp(size_t x, size_t n, const double* data);

/// @brief hill kinetics (e.g. 2X -> ... has rate n(n-1) / (2c+n(n-1)) for disassociation constant c)
/// @param[in] x the current copy number
/// @param[in] n the multiplier
/// @param[in] data the disassociation constant
double crntk_kinetics_hill(size_t x, size_t n, const double* data);




/* Lastly, we implement all of the matrix operations */

/// @brief      Computes y=Ax, where A is from the chemical master equation
/// @param[in]  crn an opaque pointer
/// @param[in]  x   the input vector
/// @param[out] y   the output vector
void crntk_id_apply(const crntk crn, const double* x, double* y);

/// @brief      Computes y=xA, where A is from the chemical master equation
/// @param[in]  crn an opaque pointer
/// @param[in]  x   the input vectpr
/// @param[out] y   the output vector
void crntk_tr_apply(const crntk crn, const double* x, double* y);

/// @brief Computes x=inv(w*L-D)*x, given tr(A)=L-D+U
/// @param[in]    crn an opaque pointer
/// @param[in]    w   the relaxation parameter
/// @param[inout] x   the input (and output) vector
void crntk_tr_sor_forward(const crntk crn, double w, double* x);

/// @brief Computes x=inv(w*U-D)*x, given tr(A)=L-D+U
/// @param[in]    crn an opaque pointer
/// @param[in]    w   the relaxation parameter
/// @param[inout] x   the input (and output) vector
void crntk_tr_sor_backward(const crntk crn, double w, double* x);

#ifdef __cplusplus
}
#endif

#endif




/// @example nullspace.c
/// @brief Find the stationary distribution using power-iteration and @p crntk_id_apply

/// @example extinction.c
/// @brief Compute the first two moments of the extinction time of a species using BiCGSTAB and @p crntk_tr_apply

/// @example birth.c
/// @brief Compute the first two moments of the birth time of a species using BiCGSTAB and @p crntk_tr_apply
