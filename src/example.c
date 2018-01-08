#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ATTN: we need to give any user-program a reference to the header file
//       containing the function and struct definitions
#include "../include/crntk.h"


// This vaguely-named function, which is fully defined after the main() 
// function, computes a single BiCGSTAB iteration using the 
// crntk_id_apply() function. Since we don't ever form the reaction rate 
// matrix, we need an iterative method.
// Of course, you are welcome to apply crntk_id_apply() to each of the 
// canonical basis vectors to explicitly form the matrix, then QR-factorize
// it, or whatever. Entirely up to you!
void solver(const crntk_crn*, const double*, double*, double*, double*, double*, double*, double*);


// ATTN: read through this main() function to get an understanding of what *you* need to do
//       to setup CRNTK to solve your problem.
int main(int argc, char** argv) {
	// call the program like ./bin/example subtrate enzyme tolerance
	size_t substrate = atoi(argv[1]);
	size_t enzyme    = atoi(argv[2]);
	double tolerance = atof(argv[3]);

	// Let's study S + E <=> X -> P + E and P -> 0 -> S
	size_t A[] = { // each column is a conservation law
		1, 0, // S (subtrate)
		1, 0, // P (product)
		1, 1, // X (dimer)
		0, 1, // E (enzyme)
		1, 0  // SLACK
	};
	size_t b[] = { // the values for each law
		substrate, enzyme
	};

	crntk_reactant reactants[] = {
		// note, these use the same ordering as the conservation laws above
		{ .affinity = &crntk_kinetics_mass_action, .data = NULL }, // the real reactants all have normal mass-action
		{ .affinity = &crntk_kinetics_mass_action, .data = NULL },
		{ .affinity = &crntk_kinetics_mass_action, .data = NULL },
		{ .affinity = &crntk_kinetics_mass_action, .data = NULL },
		{ .affinity = &crntk_kinetics_heaviside,   .data = NULL }  // assign slack heaviside dynamics to get a finite-buffer approximation
	};

	size_t complexes[] = {
	// each row is a complex (confusing, since the conservation laws were specified column-wise)
	// and of course, the ordering of reactants must be consistent with prior declarations
	//  S, P, X, E, SLACK
		1, 0, 0, 1, 0,
		0, 0, 1, 0, 0,
		0, 1, 0, 1, 0,
		0, 1, 0, 0, 0,
		0, 0, 0, 0, 1,
		1, 0, 0, 0, 0
	};

	// we have hard-coded the reaction rates here, but there is of course no need to
	// the complexes are specified by pointer offsets into the complexes[] array
	// note that complexes+5*x corresponds to the xth row, or xth complex, since there are 5 reactants
	crntk_reaction reactions[] = {
		{ .rate = 1.0,  .lhs = complexes+5*0, .rhs = complexes+5*1 }, // S+E -> X
		{ .rate = 10.0, .lhs = complexes+5*1, .rhs = complexes+5*0 }, // X -> S+E
		{ .rate = 10.0, .lhs = complexes+5*1, .rhs = complexes+5*2 }, // X -> P+E
		{ .rate = 1.0,  .lhs = complexes+5*3, .rhs = complexes+5*4 }, // P -> 0
		{ .rate = 10.0, .lhs = complexes+5*4, .rhs = complexes+5*5 }  // 0 -> S
	};

	// now we just wire everything together into a convenience object
	crntk_crn the_crn = {
		.n_reactants   = 5,
		.n_complexes   = 6,
		.n_reactions   = 5,
		.n_constraints = 2,
		.reactants = reactants,
		.complexes = complexes,
		.reactions = reactions,
		.constraints = A,
		.constraint_values = b,
	};

	// this computes the state-space and sets up the index_of functionality
	crntk_init(&the_crn);






	// now that's it! we can call crntk_id_apply(&the_crn, x, y) or whatever
	// and it will work as expected. In this example, We are going to compute 
	// the nullspace of the CME matrix to find the ergodic distribution.
	// The BiCGSTAB iterations will require some memory buffers...
	double *r0 = calloc(the_crn.n_states, sizeof(double));
	double *x = calloc(the_crn.n_states, sizeof(double));
	double *r = calloc(the_crn.n_states, sizeof(double));
	double *p = calloc(the_crn.n_states, sizeof(double));
	double *s = calloc(the_crn.n_states, sizeof(double));
	double *temp1 = calloc(the_crn.n_states, sizeof(double));
	double *temp2 = calloc(the_crn.n_states, sizeof(double));
	// along with a properly-configured initial guess
	x[0] = 1.0;
	crntk_id_apply(&the_crn, x, s);
	s[0] = 1.0;
	for (size_t i = 1; i < the_crn.n_states; i++) {
		p[i] = r0[i] = r[i] = -s[i];
	}
	
	// now we iterate until we're sufficiently close to the nullspace
	double maxnorm;
	do {
		maxnorm = 0.0;
		solver(&the_crn, r0, x, r, p, s, temp1, temp2);
		for (size_t i = 0; i < the_crn.n_states; i++) {
			maxnorm = fmax(maxnorm, fabs(r[i]));
		}
	} while (maxnorm > tolerance);

	// now we marginalize the ergodic distribution to get the enzyme copy number
	double *pmf = calloc(enzyme+1, sizeof(double));
	for (size_t i = 0; i < the_crn.n_states; i++) {
		// how much enzyme in the current state?
		size_t num_enzyme = the_crn.states[5*i+3];
		// add the probability mass to the appropriate entry
		pmf[num_enzyme] += x[i];
	}

	for (size_t i = 0; i <= enzyme; i++) {
		printf("p(%zu) = %f\n", i, pmf[i]);
	}

	// we are done, so clean up our memory buffers
	free(x); free(r); free(r0); free(p); free(s); free(temp1); free(temp2); free(pmf);

	// ATTN: be a good citizen and be sure to release your CRN.
	crntk_destroy(&the_crn);

	return 0;
}







// We wish to find the nullspace of some matrix A, where p'=Ap. We know a priori the nullspace is 1-dimensional.
// Let's replace the first row of A with that of all 1's to form B, then solve Bx=e_1.
// We use unpreconditioned bicgstab since it is fairly easy to implement and sufficiently fast
// Algorithm 7.6, https://web.stanford.edu/class/cme324/saad.pdf
void solver(const crntk_crn* crn, const double *r0, double *x, double *r, double *p, double *s, double *Ap, double *As) {
	// set temp <- A*p_id, then replace the first entry with the dot product of all 1's
	crntk_id_apply(crn, p, Ap);
	Ap[0] = 0.0;
	for (size_t i = 0; i < crn->n_states; i++) {
		Ap[0] += p[i];
	}

	// compute alpha
	double numer = 0.0, denom = 0.0;
	for (size_t i = 0; i < crn->n_states; i++) {
		numer += r0[i] * r[i];
		denom += r0[i] * Ap[i];
	}
	double alpha = numer / denom;

	// initialize s and multiply it by A
	for (size_t i = 0; i < crn->n_states; i++) {
		s[i] = r[i] - alpha * Ap[i];
	}
	crntk_id_apply(crn, s, As);
	As[0] = 0.0;
	for (size_t i = 0; i < crn->n_states; i++) {
		As[0] += s[i];
	}


	double wnumer = 0.0;
	double wdenom = 0.0;
	for (size_t i = 0; i < crn->n_states; i++) {
		wnumer += As[i] * s[i];
		wdenom += As[i] * As[i];
	}
	double omega = wnumer / wdenom;
	for (size_t i = 0; i < crn->n_states; i++) {
		x[i] += alpha * p[i] + omega * s[i];
		r[i] = s[i] - omega * As[i];
	}

	denom = numer;
	numer = 0.0;
	for (size_t i = 0; i < crn->n_states; i++) {
		numer += r[i] * r0[i];
	}

	double beta = numer / denom * alpha / omega;
	for (size_t i = 0; i < crn->n_states; i++) {
		p[i] = r[i] + beta * (p[i] - omega*Ap[i]);
	}
}
