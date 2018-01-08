# crntk
A C99-compliant implementation of stochastic chemical reaction network (chemical master equation) operations.

## Purpose
This library provides a straight-forward interface to specify a chemical reaction network and its conservation laws.
The corresponding chemical master equation, $p'=Ap$, is then exposed to the user via a few operations. In particular,
*matrix-free* evaluations of $Ax$ and $xA$ are made available, along with a few preconditioners. Because of the matrix-free
formulation, the memory-requirements are independent of the number of reactions in the network. Additionally, a novel
combinatoric formula is used to efficiently look-up array offsets for copy number configurations -- effectively a perfect 
hash map with O(1) lookup properties.

## Documentation
This library is under heavy development at the moment. For now, only the C99 interface is provided. The intent is that
in the near future, more elegant Python and Matlab wrappers will be made available. Nonetheless, the C99 interface is entirely
self-contained, with no external dependencies. While my focus is on developing these wrappers, the example file (src/example.c) 
and the header file (include/crntk.h) are both heavily commented and provide a full overview of the library and its capabilities.
Please feel free to email me with questions if you have any.

## Notes
To get running quickly, from any terminal (OSX, Linux, Windows+WSL, or Windows+msys):
> git clone https://github.com/jasondark/crntk

> cd crntk

> make example

> ./bin/example 20 5 1e-4


If that works, you're in business. The Makefile uses gcc by default and also compiles with -fopenmp:
the crntk_id_apply() and crntk_tr_apply() methods are fully parallelized. If your state-space is very small,
your program might run faster by disabling Openmp (e.g. removing that flag).

## Funding Acknowledgment
The core ideas of this software were developed while I was a graduate student at the University of California, Merced, where I received both intramural funding from the applied mathematics department and extramural funding from an NSF grant. The current implementation was developed while a postdoc at the University of California, Irvine, where my salary was derived from a DARPA grant. The content is solely the responsibility of the author and does not necessarily represent the official views of either of these agencies.
