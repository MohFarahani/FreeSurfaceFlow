#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include <math.h>

class SORSolver
{
public:

    // Constructor to manually create SORSolver
    SORSolver(int i_max, int j_max, int itermax, real ep, real omega, int norm, int checkf);

    // Constructor to create a SORSolver from a parsed configuration file
    SORSolver(const FileReader &configuration);


    // solve the pressure equation on the staggered grid
    bool solve(StaggeredGrid &grid);

    // Assignment Operator
    inline SORSolver &operator = (const SORSolver &s);

private:
    // solver values
    int imax, jmax, smax, normfreq, checkfreq;
    real eps, om;
};



#endif //SOR_SOLVER_HH

