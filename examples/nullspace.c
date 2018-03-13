#include <stdio.h>
#include <stdlib.h>
#include <math.h> // fabs, fmax

#include "../include/crntk.h"


int main(int argc, char** argv) {
    if (argc < 5) {
        printf("Usage: ./nullspace mass kf kb error\n");
        return 0;
    }

    // get the arguments
    const size_t n = atoi(argv[1]);
    const double kf = atof(argv[2]);
    const double kb = atof(argv[3]);
    const double tol = atof(argv[4]);

    // Consider the chemical reaction network 2A <=> B (2A->B with rate kf, B->2A with rate kb)
    // Intuition demands a unique stationary distribution for the joint copy numbers.
    // We specify the chemical reaction network below then compute this distribution.

    // Initialize the system
    crntk demo;
    if (!crntk_init(&demo, 2, 2, 2, 1)) {
        printf("Could not allocate memory for crntk stage 1\n");
        return 0;
    }

    // reactants: the ordering here persists across all other API calls
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // A
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // B

    size_t cplx[2];

    // constraints: just one conservation law at work here
    cplx[0] = 1; cplx[1] = 2;
    crntk_add_constraint(demo, n, cplx); // n = A + 2B

    // complexes: just the two
    cplx[0] = 2; cplx[1] = 0;
    crntk_add_complex(demo, cplx); // 2*A
    cplx[0] = 0; cplx[1] = 1;
    crntk_add_complex(demo, cplx); // B

    // reactions: again, just the one reversible reaction (represented as two separate reactions)
    crntk_add_reaction(demo, kf, 0, 1); // 2A -> B
    crntk_add_reaction(demo, kb, 1, 0); // B -> 2A

    // We need to finalize the CRN pointer before we can do anything with it.
    // Warning: Depending on the constraints, this might allocate A LOT of memory.
    if (!crntk_finalize(demo)) {
        printf("Could not allocate memory for crntk stage 2\n");
        return 0;
    }


    // The rest of this code is just numerical linear algebra, completely independent
    // of the CRN stuff above. Go to line 135 to see how the program output is generated.

    // Goal: Find a vector in the nullspace of A, where the CME reads p'=Ap.
    // Observations: All eigenvalues of A are non-positive by construction.
    // Calculation: 1. Apply Gershgorin's circle theorem to shift matrix by a positive amount.
    //                 This makes the least dominant eigenvector (the zero one) the most dominant.
    //              2. Find this eigenvector by power iteration.

    // we initialize some memory to store the intermediate vectors
    const size_t len = crntk_dim(demo); // the number of unique states
    double *x0 = calloc(len, sizeof(double)); // two intermediate vectors
    double *x1 = calloc(len, sizeof(double));
    // scratch variables
    double *swap;
    double temp;
    double norm;

    // get the diagonal elements and find the biggest one (most negative)
    double shift = 0.0;
    crntk_diag(demo, x1);
    for (size_t i = 0; i < len; i++) {
        shift = fmin(shift, x1[i]);
    }
    shift *= -2.0;

    // Power iteration: shift the matrix, find the new dominant eigenvector
    x0[0] = 1.0; // initial guess: x0=(1,0,0,...)
    do  {
        // Compute x1 = Ax0 (note, we really want x1 = (A+shift*I)x0)
        crntk_id_apply(demo, x0, x1);
        // Here we correct the +shift*x0 bit, as well as compute the L1 normalization
        temp = 0.0;
        for (size_t i = 0; i < len; i++) {
            x1[i] += shift * x0[i];
            temp += fabs(x1[i]);
        }
        // normalize and compute maxnorm for convergence check
        temp = 1.0 / temp;
        x1[0] *= temp;
        norm = fabs(x0[0] - x1[0]);
        for (size_t i = 1; i < len; i++) {
            x1[i] *= temp;
            norm = fmax(norm, fabs(x0[i]-x1[i]));
        }
        // swap x0 and x1 for next iteration
        swap = x1;
        x1 = x0;
        x0 = swap;
    } while (norm > tol);


    // print out the joint distribution
    printf("A\tB\tp\n");
    for (size_t i = 0; i < len; i++) {
        const size_t* state = crntk_state(demo, i);
        printf("%zu\t%zu\t%f\n", state[0], state[1], x1[i]);
    }

    // we are done, so clean up our memory buffers
    free(x0); free(x1);

    // be a good citizen and be sure to release the CRN memory as well
    crntk_destroy(demo);

    return 0;
}
