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
    if (argc < 7) {
        printf("Usage: ./birth S0 E0 kf kb kcat error\n");
        return 0;
    }

    // get the arguments
    const size_t S0 = atoi(argv[1]);
    const size_t E0 = atoi(argv[2]);
    const double kf = atof(argv[3]);
    const double kb = atof(argv[4]);
    const double kcat = atof(argv[5]);
    const double tol = atof(argv[6]);

    const double km = (kb+kcat) / kf;
    // Consider the Michaelis-Menten system S+E <=> ES -> P+E
    // Perhaps we are interested in the time to the formation of the first P.
    // As we saw in extinction.c, non-physical reactions leak probability flux
    // which can be used to construct first-passage time moments.
    // An alternative to a non-physical rate (allowing a reaction to fire even if copy numbers are not large enough)
    // is a non-conservative reaction, which physically removes mass from the system.
    // That is, consider instead S+E <=> ES -> 0.

    // For fun, we will also compare this to the classic Michaelis-Menten approximation

    // Initialize the system
    crntk demo;
    if (!crntk_init(&demo, 3, 3, 3, 2)) {
        printf("Could not allocate memory for crntk stage 1\n");
        return 0;
    }

    // reactants: the ordering here persists across all other API calls
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // S
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // E
    crntk_add_reactant(demo, &crntk_kinetics_mass_action, NULL); // ES

    // constraints:
    size_t cplx[3];

    cplx[0] = 1; cplx[1] = 0; cplx[2] = 1;
    crntk_add_constraint(demo, S0, cplx); // S0 = S+ES
    cplx[0] = 0; cplx[1] = 1; cplx[2] = 1;
    crntk_add_constraint(demo, E0, cplx); // E0 = E+ES

    // complexes:
    cplx[0] = 1; cplx[1] = 1; cplx[2] = 0;
    crntk_add_complex(demo, cplx); // S+E
    cplx[0] = 0; cplx[1] = 0; cplx[2] = 1;
    crntk_add_complex(demo, cplx); // ES
    cplx[0] = 0; cplx[1] = 0; cplx[2] = 0;
    crntk_add_complex(demo, cplx); // 0

    // reactions:
    crntk_add_reaction(demo, kf, 0, 1); // Complex 0 -> Complex 1
    crntk_add_reaction(demo, kb, 1, 0); // Complex 1 -> Complex 0
    crntk_add_reaction(demo, kcat, 1, 2); // Complex 1 -> Complex 2

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
    const size_t init[] = { S0, E0, 0 }; // (S,E,ES)
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

    printf("Exact:\n");
    printf("Birth time (mean):\t%f\n", mu1);
    printf("Birth time (stdv):\t%f\n", sqrt(mu2 - mu1*mu1));

    // We know that if S converted at some rate implies the formation time
    // is exponentially distributed with this rate.
    double rate = (kcat * (double) E0 * (double) S0) / (km + (double) S0);
    printf("\nApproximate: Exp(lambda=%f)\n", rate);
    printf("Birth time (mean):\t%f\n", 1.0/rate);
    printf("Birth time (stdv):\t%f\n", 1.0/rate);


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
