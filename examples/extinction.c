#include <stdio.h>
#include <stdlib.h>
#include <math.h> // fabs, fmax

#include "../include/crntk.h"

// declare the bicgstab solver
size_t bicgstab(const void*, void (*)(const void*, const double*, double*), const size_t, const double *, double *, double, size_t);


// an implementation of matrix-multiplication with a preconditioner
void matmul(const void* data, const double *x, double *y) {
    const crntk crn = (crntk) data;

    // apply tr(A)
    crntk_tr_apply(crn, x, y);

    // then the left preconditioner. choices are:
    // crntk_tr_sor_backward, crntk_tr_sor_forward, crntk_tr_ssor_backward, crntk_tr_ssor_forward
    // if no preconditioner is used, we need to negate the vector y.
    crntk_tr_sor_forward(crn, 1.0, y); // 1 forward Jacobi iteration as a preconditioner
    // note that at the moment, none of the preconditioners are parallelized
}



int main(int argc, char** argv) {
    if (argc < 5) {
        printf("Usage: ./extinction mass a b error\n");
        return 0;
    }

    // get the arguments
    const size_t n = atoi(argv[1]);
    const double a = atof(argv[2]);
    const double b = atof(argv[3]);
    const double tol = atof(argv[4]);

    // Consider the chemical reaction network X+A -> 2X (with rate a) and
    //                                          X -> A (with rate b).
    // As the total conserved mass X+A -> infinity, the extinction event X=0 happens
    // after increasingly long periods of time. Let's explore this extinction time.

    // To do so, we introduce a "slack reactant" to enforce A <= n-1. That is, consider a Z such that
    // A+Z = n-1. We modify the original reactions as follows:
    // X+A -> 2X+Z (rate a) and X+Z -> A (rate b)
    // We also endow Z with non-physical, constant kinetics so as not to change the reaction propensities.
    // This permits probability flux to "leak" out of the state (X,A,Z) = (1,n-1,0).
    // The mean extinction time is then <p0, -inv(tr(A))*1>.
    // (See Chou, T. and M. R. D'Orsogna (2014).
    //  "First-passage phenomena and their applications".
    //   In: ed. by Ralf Metzler, Gleb Oshanin, and Sidney Redner.
    //   Vol. 35. World Scientific. Chap. 13, pp. 306–345.)


    // Initialize the system
    crntk demo;
    if (!crntk_init(&demo, 3, 4, 2, 2)) {
        printf("Could not allocate memory for crntk stage 1\n");
        return 0;
    }

    // reactants: the ordering here persists across all other API calls
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // X
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // A
    crntk_add_reactant(demo, &crntk_kinetics_fsp, NULL);         // Z

    // constraints: just one conservation law at work here
    crntk_add_constraint(demo, n,   1, 1, 0); //   n = X+A
    crntk_add_constraint(demo, n-1, 0, 1, 1); // n-1 = A+Z

    // complexes: each side of each reaction
    crntk_add_complex(demo, 1, 1, 0); // X+A
    crntk_add_complex(demo, 2, 0, 1); // 2X+Z
    crntk_add_complex(demo, 1, 0, 1); // X+Z
    crntk_add_complex(demo, 0, 1, 0); // A

    // reactions: again, just the one reversible reaction (represented as two separate reactions)
    crntk_add_reaction(demo, a, 0, 1); // Complex 0 -> Complex 1 (X+A -> 2X+Z)
    crntk_add_reaction(demo, b, 2, 3); // Complex 2 -> Complex 3 (X+Z -> A)

    // We need to finalize the CRN pointer before we can do anything with it.
    // Warning: Depending on the constraints, this might allocate A LOT of memory.
    if (!crntk_finalize(demo)) {
        printf("Could not allocate memory for crntk stage 2\n");
        return 0;
    }


    // we initialize some memory to store the intermediate vectors
    const size_t len = crntk_dim(demo); // the number of unique states
    double *rhs = calloc(len, sizeof(double));
    double *x = calloc(len, sizeof(double));

    // suppose we want only the extinction time statistics for a particular state, e.g.
    const size_t init[] = { n, 0, n-1 }; // X=n, A=0, Z=n-1
    const size_t index = crntk_index_of(demo, init);
    // x[index] will have the statistics of interest


    // Setup solve for first moment
    for (size_t i = 0; i < len; i++) {
        rhs[i] = -1.0;
    }
    crntk_tr_sor_forward(demo, 1.0, rhs); // important that we use the same preconditioner here as above
    bicgstab((void*) demo, matmul, len, rhs, x, tol, 10000);

    double mu1 = x[index];
    
    // Setup solve for second (uncentered) moment
    for (size_t i = 0; i < len; i++) {
        rhs[i] = -2.0 * x[i];
        x[i] = 0.0;
    }
    crntk_tr_sor_forward(demo, 1.0, rhs); // important that we use the same preconditioner here as above
    bicgstab((void*) demo, matmul, len, rhs, x, tol, 10000);

    double mu2 = x[index];


    printf("Extinction time (mean):\t%f\n", mu1);
    printf("Extinction time (stdv):\t%f\n", sqrt(mu2 - mu1*mu1));

    // we are done, so clean up our memory buffers
    free(x); free(rhs);

    // be a good citizen and be sure to release the CRN memory as well
    crntk_destroy(demo);

    return 0;
}




// We need some iterative scheme to solve Ax=b. BiCGSTAB is easy enough to implement, so here it is.
// Algorithm 7.6, https://web.stanford.edu/class/cme324/saad.pdf
size_t bicgstab(
    const void* data,
    void (*apply)(const void*, const double*, double*),
    const size_t n,
    const double *b,
    double *x,
    double tol,
    size_t maxiter)
{
    size_t niter = maxiter;

    double *r0 = calloc(n, sizeof(double));
    double *r = calloc(n, sizeof(double));
    double *p = calloc(n, sizeof(double));
    double *s = calloc(n, sizeof(double));
    double *Ap = calloc(n, sizeof(double));
    double *As = calloc(n, sizeof(double));

    apply(data, x, r0);
    for (size_t i = 0; i < n; i++) {
        p[i] = r0[i] = r[i] = b[i] - r0[i];
    }


    for (size_t iter = 0; iter < maxiter; iter++) {
        apply(data, p, Ap);

        // compute alpha
        double numer = 0.0, denom = 0.0;
        for (size_t i = 0; i < n; i++) {
            numer += r0[i] * r[i];
            denom += r0[i] * Ap[i];
        }
        double alpha = numer / denom;

        // initialize s and multiply it by A
        for (size_t i = 0; i < n; i++) {
            s[i] = r[i] - alpha * Ap[i];
        }
        apply(data, s, As);

        double wnumer = 0.0;
        double wdenom = 0.0;
        for (size_t i = 0; i < n; i++) {
            wnumer += As[i] * s[i];
            wdenom += As[i] * As[i];
        }
        double omega = wnumer / wdenom;

        double maxnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * As[i];
            maxnorm = fmax(maxnorm, fabs(r[i]));
        }
        if (maxnorm < tol) {
            niter = iter;
            break;
        }
    
        denom = numer;
        numer = 0.0;
        for (size_t i = 0; i < n; i++) {
            numer += r[i] * r0[i];
        }

        double beta = numer / denom * alpha / omega;
        for (size_t i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega*Ap[i]);
        }
    }

    free(r0); free(r); free(p); free(s); free(Ap); free(As);
    return niter;
}
