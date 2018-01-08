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

#include "../include/crntk.h"

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
/*double crntk_kinetics_hill(size_t x, size_t n, const double *data) {
    if (x < n) {
        return 0.0;
    }
    double kd = data[0];
    double coef = data[1];
    double z = 1.0;
    for (size_t i = 0; i < n; i++) {
        double d = pow((double) (x-i), coef);
        z *= d / (kd + d);
    }
    return z;
}*/


static inline bool react(const crntk_crn *crn, const size_t rxn, const size_t *n0, size_t *n1, double *propensity) {
    const size_t *lhs = crn->reactions[rxn].lhs;
    const size_t *rhs = crn->reactions[rxn].rhs;

    bool valid = true;
    *propensity = crn->reactions[rxn].rate;

    // iterate through each species
    for (size_t k = 0; k < crn->n_reactants; k++) {
        // check for sufficient copy number
        valid = valid && (n0[k] >= lhs[k]);
        *propensity *= crn->reactants[k].affinity(n0[k], lhs[k], crn->reactants[k].data);
    }

    // update the n1 state (if valid==false, this is nonsense)
    for (size_t k = 0; k < crn->n_reactants; k++) {
        n1[k] = n0[k] + rhs[k] - lhs[k];
    }

    return valid;
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
static bool init_lattice_table(crntk_crn *lattice) {
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
static bool init_lattice_states(crntk_crn *lattice) {
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

size_t crntk_index_of(const crntk_crn *lattice, const size_t *index) {
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
// The next 3 functions are for initializing the table field of a lattice struct





// y = id(A)x
void crntk_id_apply(const crntk_crn* crn, const double *x, double *y) {
    // init y to zero vector
    memset(y, 0, crn->n_states * sizeof(double));

    #pragma omp parallel for
    for (size_t i = 0; i < crn->n_states; i++) {
        // n0 is the current state, n1 is a temporary variable for neighboring states
        const size_t *n0 = crn->states + crn->n_reactants * i;
        size_t n1[crn->n_reactants];

        // a temporary variable for the reaction propensity and neighbors
        double propensity;
        bool nbrd;

        for (size_t k = 0; k < crn->n_reactions; k++) {
            nbrd = react(crn, k, n0, n1, &propensity);
            propensity *= x[i];
            y[i] -= propensity;
            if (nbrd) {
                const size_t j = crntk_index_of(crn, n1);
                #pragma omp atomic
                y[j] += propensity;
            }
        }
    }
}



// y = tr(A)x
void crntk_tr_apply(const crntk_crn* crn, const double *x, double *y) {
    // init y to zero vector
    memset(y, 0, crn->n_states * sizeof(double));

    #pragma omp parallel for
    for (size_t i = 0; i < crn->n_states; i++) {
        // n0 is the current state, n1 is a temporary variable for neighboring states
        const size_t *n0 = crn->states + crn->n_reactants * i;
        size_t n1[crn->n_reactants];

        // a temporary variable for the reaction propensity and neighbors
        double propensity;
        bool nbrd;

        for (size_t k = 0; k < crn->n_reactions; k++) {
            nbrd = react(crn, k, n0, n1, &propensity);
            y[i] += propensity * ((nbrd? x[crntk_index_of(crn, n1)] : 0.0) - x[i]);
        }
    }
}


// Given, tr(A) = L+U-D...
// (1) compute x <- inv(w*L-D) x
void crntk_tr_sor_forward(const crntk_crn* crn, double omega, double *x) {
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
            nbrd = react(crn, k, n0, n1, &propensity);
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
void crntk_tr_sor_backward(const crntk_crn* crn, double omega, double *x) {
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
            nbrd = react(crn, k, n0, n1, &propensity);
            flux -= propensity;
            j = nbrd? crntk_index_of(crn, n1) : i-1;
            if (j > i) {
                x[i] -= omega * propensity * x[j];
            }
        }
        x[i] /= flux;
    }
}

static inline void jacobi(const crntk_crn* crn, double *x) {
    #pragma omp parallel for
    for (size_t i = 0; i < crn->n_states; i++) {
        // n0 is the current state, n1 is a temporary variable for neighboring states
        const size_t *n0 = crn->states + crn->n_reactants * i;
        size_t n1[crn->n_reactants];

        // a temporary variable for the reaction propensity and neighbors
        double flux = 0.0;
        double propensity;

        for (size_t k = 0; k < crn->n_reactions; k++) {
            react(crn, k, n0, n1, &propensity);
            flux -= propensity;
        }
        x[i] *= flux;
    }
}


void crntk_tr_ssor_forward(const crntk_crn* crn, double omega, double *x) {
    crntk_tr_sor_forward(crn, omega, x);
    jacobi(crn, x);
    crntk_tr_sor_backward(crn, omega, x);
}
void crntk_tr_ssor_backward(const crntk_crn* crn, double omega, double *x) {
    crntk_tr_sor_backward(crn, omega, x);
    jacobi(crn, x);
    crntk_tr_sor_forward(crn, omega, x);
}







// and now for the lattice stuff


// version 1 (now): brute force enumerate the states (with some clever shortcuts)
// version 2 (later): it is possible to construct these efficiently on the fly
//            by reducing A to column Hermite Normal form...
//            This is a possible memory/speed trade-off to be considered later.

bool crntk_init(crntk_crn *lattice) {
	if (!init_lattice_table(lattice))
		return false;

	if (!init_lattice_states(lattice))
		return false;

	return true;
}
void crntk_destroy(crntk_crn *lattice) {
	for (size_t i = 0; i < lattice->n_reactants; i++)
		free(lattice->table[i]);
	free(lattice->table);
	free(lattice->states);
}
