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

#include <stdarg.h> // varargs for add_complex and add_constraint
#include <stdlib.h> // calloc
#include <string.h> // memset, memcpy


#include "../include/crntk.h"

// internal struct: a reactant object
typedef struct {
    double (*affinity)(size_t,size_t,const double*);
    double *data;
} crntk_reactant;

// internal struct: a reaction object
typedef struct {
    double rate; // the rate constant for the reaction
    size_t* lhs; // a pointer to the source complex
    size_t* rhs; // a pointer to the destination complex
    bool conservative; // does the reaction satisfy the conservation laws?
} crntk_reaction;

/// implementation of crntk struct
struct crntk_crn {
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

    // internal builder quantities
    crntk_reactant* _ptr_reactant;
    size_t* _ptr_complex;
    crntk_reaction* _ptr_reaction;
    size_t* _ptr_constraint;
    size_t* _ptr_value;
};


// the next 5 functions are internal helper functions
static inline void react(
    // input parameters
    const crntk *crn, const size_t rxn_index, const size_t *n0,
    // output parameters
    size_t *n1, double *propensity, bool *bounded)
{
    // copy the reaction
    crntk_reaction rxn = crn->reactions[rxn_index];

    // iterate through each species
    for (size_t k = 0; k < crn->n_reactants; k++) {
        // check for sufficient copy number
        rxn.conservative = rxn.conservative && (n0[k] >= rxn.lhs[k]);
        // compute the propensity regardless -- this permits non-physical kinetics like FSP or non-conservative reactions
        rxn.rate *= crn->reactants[k].affinity(n0[k], rxn.lhs[k], crn->reactants[k].data);
    }

    // update the output parameters
    for (size_t k = 0; k < crn->n_reactants; k++) {
        n1[k] = n0[k] + rxn.rhs[k] - rxn.lhs[k];
    }
    *propensity = rxn.rate;
    *bounded = rxn.conservative;
}
static inline size_t tensor_offset(size_t *index, const size_t rank, size_t *bound) {
    size_t offset = 0;
    for (size_t i = 0; i < rank; i++) {
        offset = *index + offset * (1 + *bound);
        bound += 1;
        index += 1;
    }
    return offset;
}
static inline size_t tensor_delta(size_t *b, size_t *c, size_t *n, const size_t rank) {
    size_t i = 0, min, temp;
    for (; i < rank && n[i] == 0; i++);
    min = (b[i] - c[i]) / n[i];
    for (i++; i < rank; i++) {
        if (n[i] == 0) continue;
        temp = (b[i] - c[i]) / n[i];
        if (temp < min) min = temp;
    }
    return min;
}
static bool init_lattice_table(crntk *lattice) {
    const size_t m = lattice->n_constraints, n = lattice->n_reactants;
    size_t *A = lattice->constraints, *b = lattice->constraint_values;

    size_t** table = calloc(n, sizeof(size_t*));
    if (table == NULL)
        return false;

    size_t numel = 1;
    for (size_t i = 0; i < m; i++)
        numel *= b[i]+1;

    for (size_t i = 0; i < n; i++) {
        table[i] = calloc(numel, sizeof(size_t));
        if (table[i] == NULL) {
            for (size_t j = 0; j < i; j++)
                free(table[j]);
            free(table);
            return false;
        }
    }

    size_t index[m];
    for (size_t i = 0; i < m; i++)
        index[i] = 0;


    // start with last variable
    A += m*(n-1);
    size_t delta = tensor_offset(A, m, b);
    size_t inc = tensor_delta(b, index, A, m);
    for (size_t j = 0; j <= inc; j++)
        table[0][j*delta] = 1;

    // and then the remaining variables
    for (size_t k = 1; k < n; k++) {
        A -= m;
        delta = tensor_offset(A, m, b);

        memcpy(index, b, m * sizeof(size_t));

        for (size_t i = 0; i < numel; i++) {
            for (size_t j = m; j-- > 0; ) {
                if (index[j] < b[j]) {
                    index[j] += 1;
                    break;
                }
                index[j] = 0;
            }

            inc = tensor_delta(b, index, A, m);

            for (size_t j = 0; j <= inc; j++)
                table[k][i+j*delta] += table[k-1][i];
        }
    }

    lattice->table = table;
    return true;
}
static bool init_lattice_states(crntk *lattice) {
    const size_t m = lattice->n_constraints, n = lattice->n_reactants;

    const size_t dim = lattice->table[n-1][tensor_offset(lattice->constraint_values, m, lattice->constraint_values)];
    lattice->n_states = dim;
    lattice->states = calloc(n * dim, sizeof(size_t));
    if (lattice->states == NULL)
        return false;

    size_t b[m], x[n];
    memcpy(b, lattice->constraint_values, m*sizeof(size_t));
    for (size_t i = 0; i < n; i++)
        x[i] = 0;

    size_t i = 0;
    size_t j = n-1;
    size_t k;
    bool flag, finished = false;

    while (!finished) {
        // Special case 1: first column
        if (j == 0) {
            for (k = 0; k < m; k++) {
                // if no more increments are possible we are done
                if (lattice->constraints[k+j*m] > b[k]) {
                    finished = true;
                    break;
                }
            }
            if (!finished) {
                x[j] += 1;
                for (k = 0; k < m; k++)
                    b[k] -= lattice->constraints[k+j*m];
                j = n-1;
            }
            continue;
        }

        // we can easily check if there are solutions. if not, let's not waste time
        if (lattice->table[n-1-j][tensor_offset(b, m, lattice->constraint_values)] == 0) {
            j -= 1;
            continue;
        }

        // Special case 2: last column
        if (j == n-1) {
            // we know there is a solution, so we record it
            for (k = 0; k < m; k++) {
                if (b[k] == 0 && lattice->constraints[k+j*m] == 0)
                    continue;
                else
                    break;
            }
            x[j] = b[k] / lattice->constraints[k+j*m];
            memcpy(lattice->states + n*i, x, n * sizeof(size_t));
            i += 1;
            x[j] = 0;
            j -= 1;
            continue;
        }



        // general case: if we can increment, great, otherwise we back-track 
        flag = true;
        for (k = 0; k < m; k++) {
            if (lattice->constraints[k+j*m] > b[k]) {
                flag = false;
                break;
            }
        }
        if (flag) {
            x[j] += 1;
            for (k = 0; k < m; k++)
                b[k] -= lattice->constraints[k+j*m];
            j = n-1;
        }
        else {
            if (x[j] != 0) {
                for (k = 0; k < m; k++)
                    b[k] += x[j] * lattice->constraints[k+j*m];
                x[j] = 0;
            }
            j -= 1;
        }
    }

    return true;
}







double crntk_kinetics_mass_action(size_t x, size_t n, const double *data) {
    if (x < n) {
        return 0.0;
    }

    double response = 1.0;
    for (size_t i = 0; i < n; i++) {
        response *= (double) (x-i) / (double) (i+1);
    }
    return response;
}
double crntk_kinetics_heaviside(size_t x, size_t n, const double *data) {
    return (x < n)? 0.0 : 1.0;
}
double crntk_kinetics_fsp(size_t x, size_t n, const double *data) {
    return 1.0;
}
double crntk_kinetics_hill(size_t x, size_t n, const double *data) {
    const double temp = crntk_kinetics_mass_action(x, n, NULL);
    return temp / (data[0] + temp);
}



size_t crntk_index_of(const crntk *lattice, const size_t *index) {
    const size_t m = lattice->n_constraints, n = lattice->n_reactants;
    const size_t nullity = n-m;
    size_t *A = lattice->constraints;

    size_t offset = 0;
    
    size_t b[m];
    memcpy(b, lattice->constraint_values, m*sizeof(size_t));

    for (size_t i = 0; i < nullity; i++) {
        for (size_t j = 0; j < index[i]; j++) {
            offset += lattice->table[n-i-2][tensor_offset(b, m, lattice->constraint_values)];
            for (size_t k = 0; k < m; k++)
                b[k] -= A[k];
        }
        A += m;
    }

    return offset;
}
size_t crntk_dim(const crntk *crn) {
    return crn->n_states;
}

const size_t* crntk_state(const crntk* crn, const size_t i) {
    return crn->states + crn->n_reactants*i;
}

void crntk_diag(const crntk* crn, double* d) {
    #pragma omp parallel
    {
        size_t n1[crn->n_reactants];
        double propensity;
        bool valid;
    
        #pragma omp for
        for (size_t i = 0; i < crn->n_states; i++) {
            // n0 is the current state, n1 is a temporary variable for neighboring states
            const size_t *n0 = crn->states + crn->n_reactants * i;

            // a temporary variable for the reaction propensity and neighbors
            double flux = 0.0;

            for (size_t k = 0; k < crn->n_reactions; k++) {
                react(crn, k, n0, n1, &propensity, &valid);
                flux -= propensity;
            }
            d[i] = flux;
        }
    }
}


// y = id(A)x
void crntk_id_apply(const crntk* crn, const double *x, double *y) {
    // init y to zero vector
    memset(y, 0, crn->n_states * sizeof(double));

    #pragma omp parallel
    {
        // thread-local scratch variables
        size_t n1[crn->n_reactants];
        double propensity;
        bool valid;

        #pragma omp for
        for (size_t i = 0; i < crn->n_states; i++) {
            // n0 is the current state, n1 is a temporary variable for neighboring states
            const size_t *n0 = crn->states + crn->n_reactants * i;

            for (size_t k = 0; k < crn->n_reactions; k++) {
                react(crn, k, n0, n1, &propensity, &valid);
                propensity *= x[i];
                #pragma omp atomic
                y[i] -= propensity;
                if (valid) {
                    const size_t j = crntk_index_of(crn, n1);
                    #pragma omp atomic
                    y[j] += propensity;
                }
            }
        }
    }
}



// y = tr(A)x
void crntk_tr_apply(const crntk* crn, const double *x, double *y) {
    // init y to zero vector
    memset(y, 0, crn->n_states * sizeof(double));

    #pragma omp parallel
    {
        // thread-local scratch variables
        size_t n1[crn->n_reactants];
        double propensity;
        bool valid;

        #pragma omp for
        for (size_t i = 0; i < crn->n_states; i++) {
            // n0 is the current state, n1 is a temporary variable for neighboring states
            const size_t *n0 = crn->states + crn->n_reactants * i;

            for (size_t k = 0; k < crn->n_reactions; k++) {
                react(crn, k, n0, n1, &propensity, &valid);
                y[i] += propensity * ((valid? x[crntk_index_of(crn, n1)] : 0.0) - x[i]);
            }
        }
    }
}


// Given, tr(A) = L+U-D...
// (1) compute x <- inv(w*L-D) x
void crntk_tr_sor_forward(const crntk* crn, double omega, double *x) {
    size_t n1[crn->n_reactants];
    for (size_t i = 0; i < crn->n_states; i++) {
        // n0 is the current state, n1 is a temporary variable for neighboring states
        const size_t *n0 = crn->states + crn->n_reactants * i;

        // a temporary variable for the reaction propensity and neighbors
        double flux = 0.0;
        double propensity;
        bool nbrd;
        size_t j;

        for (size_t k = 0; k < crn->n_reactions; k++) {
            react(crn, k, n0, n1, &propensity, &nbrd);
            flux -= propensity;
            j = nbrd? crntk_index_of(crn, n1) : i+1;
            if (j < i) {
                x[i] -= omega * propensity * x[j];
            }
        }
        x[i] /= flux;
    }
}
// (2) compute x <- inv(w*U-D) x
void crntk_tr_sor_backward(const crntk* crn, double omega, double *x) {
    size_t n1[crn->n_reactants];
    for (size_t i = crn->n_states; i-- > 0;) {
        // n0 is the current state, n1 is a temporary variable for neighboring states
        const size_t *n0 = crn->states + crn->n_reactants * i;

        // a temporary variable for the reaction propensity and neighbors
        double flux = 0.0;
        double propensity;
        bool nbrd;
        size_t j;

        for (size_t k = 0; k < crn->n_reactions; k++) {
            react(crn, k, n0, n1, &propensity, &nbrd);
            flux -= propensity;
            j = nbrd? crntk_index_of(crn, n1) : i-1;
            if (j > i) {
                x[i] -= omega * propensity * x[j];
            }
        }
        x[i] /= flux;
    }
}

static inline void jacobi(const crntk* crn, double *x) {
    #pragma omp parallel
    {
        size_t n1[crn->n_reactants];
        double propensity;
        bool valid;
    
        #pragma omp for
        for (size_t i = 0; i < crn->n_states; i++) {
            // n0 is the current state, n1 is a temporary variable for neighboring states
            const size_t *n0 = crn->states + crn->n_reactants * i;

            // a temporary variable for the reaction propensity and neighbors
            double flux = 0.0;

            for (size_t k = 0; k < crn->n_reactions; k++) {
                react(crn, k, n0, n1, &propensity, &valid);
                flux -= propensity;
            }
            x[i] *= flux;
        }
    }
}







// and now for the lattice stuff


// version 1 (now): brute force enumerate the states (with some clever shortcuts)
// version 2 (later): it is possible to construct these efficiently on the fly
//            by reducing A to column Hermite Normal form...
//            This is a possible memory/speed trade-off to be considered later.


bool crntk_init(crntk **out, const size_t n_reactants, const size_t n_complexes, const size_t n_reactions, const size_t n_constraints) {
    crntk *crn = calloc(1, sizeof(crntk));
    crn->n_reactants   = n_reactants;
    crn->n_complexes   = n_complexes;
    crn->n_reactions   = n_reactions;
    crn->n_constraints = n_constraints;

    crn->reactants = calloc(n_reactants, sizeof(crntk_reactant));
    if (crn->reactants == NULL) {
        return false;
    }
    crn->complexes = calloc(n_complexes*n_reactants, sizeof(size_t));
    if (crn->complexes == NULL) {
        free(crn->reactants);
        return false;
    }
    crn->reactions = calloc(n_reactions, sizeof(crntk_reaction));
    if (crn->reactions == NULL) {
        free(crn->reactants);
        free(crn->complexes);
        return false;
    }
    crn->constraints = calloc(n_constraints*n_reactants, sizeof(size_t));
    if (crn->constraints == NULL) {
        free(crn->reactants);
        free(crn->complexes);
        free(crn->reactions);
        return false;

    }
    crn->constraint_values = calloc(n_constraints, sizeof(size_t));
    if (crn->constraint_values == NULL) {
        free(crn->reactants);
        free(crn->complexes);
        free(crn->reactions);
        free(crn->constraints);
        return false;
    }

    crn->_ptr_reactant   = crn->reactants;
    crn->_ptr_complex    = crn->complexes;
    crn->_ptr_reaction   = crn->reactions;
    crn->_ptr_constraint = crn->constraints;
    crn->_ptr_value      = crn->constraint_values;

    *out = crn;

    return true;
}

void crntk_add_reactant(crntk *crn, double (*affinity)(size_t,size_t,const double*), double *data) {
    crn->_ptr_reactant->affinity = affinity;
    crn->_ptr_reactant->data = data;
    crn->_ptr_reactant++;
}

void crntk_add_complex(crntk *crn, ...) {
    va_list va;
    va_start(va, crn);
    for (size_t i = 0; i < crn->n_reactants; i++) {
        *crn->_ptr_complex = va_arg(va, size_t);
        crn->_ptr_complex++;
    }
    va_end(va);
}

void crntk_add_reaction(crntk *crn, double rate, size_t lhs, size_t rhs) {
    crn->_ptr_reaction->rate = rate;
    crn->_ptr_reaction->lhs = crn->complexes + crn->n_reactants * lhs;
    crn->_ptr_reaction->rhs = crn->complexes + crn->n_reactants * rhs;
    crn->_ptr_reaction++;
}

void crntk_add_constraint(crntk *crn, size_t value, ...) {
    va_list va;
    va_start(va, value);

    *crn->_ptr_value = value;
    crn->_ptr_value++;

    for (size_t i = 0; i < crn->n_reactants; i++) {
        *crn->_ptr_constraint = va_arg(va, size_t);
        crn->_ptr_constraint += crn->n_constraints;
    }
    crn->_ptr_constraint -= crn->n_constraints * crn->n_reactants - 1;

    va_end(va);
}

bool crntk_finalize(crntk *crn) {
    // check the conservation of each reaction
    for (size_t i = 0; i < crn->n_reactions; i++) {
        crn->reactions[i].conservative = true;
        const size_t *lhs = crn->reactions[i].lhs;
        const size_t *rhs = crn->reactions[i].rhs;
        const size_t *constraint = crn->constraints;

        size_t value_lhs, value_rhs;
        for (size_t j = 0; j < crn->n_constraints; j++) {
            value_lhs = 0;
            value_rhs = 0;
            for (size_t k = 0; k < crn->n_reactants; k++) {
                value_lhs += lhs[k] * constraint[k * crn->n_constraints];
                value_rhs += rhs[k] * constraint[k * crn->n_constraints];
            }
            constraint++;

            if (value_lhs != value_rhs) {
                crn->reactions[i].conservative = false;
                break;
            }
        }
    }


    if (!init_lattice_table(crn))
        return false;

    if (!init_lattice_states(crn))
        return false;

    return true;
}

void crntk_destroy(crntk *crn) {
    free(crn->reactants);
    free(crn->complexes);
    free(crn->reactions);
    free(crn->constraints);
    free(crn->constraint_values);
    for (size_t i = 0; i < crn->n_reactants; i++)
        free(crn->table[i]);
    free(crn->table);
    free(crn->states);
    free(crn);
}
