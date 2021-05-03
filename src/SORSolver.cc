
#include "SORSolver.hh"

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

// Constructor to manually create SORSolver
SORSolver::SORSolver(int i_max, int j_max, int itermax, real ep, real omega, int norm, int checkf)
    : normfreq(1), checkfreq(1)
{
    CHECK_MSG((i_max >= 0), "wrong input for imax: " << i_max);
    CHECK_MSG((j_max >= 0), "wrong input for jmax: " << j_max);
    CHECK_MSG((itermax >= 0), "wrong input for itermax: " << itermax);
    CHECK_MSG((ep >= 0), "wrong input for eps: " << ep);
    CHECK_MSG((omega == 1 || (omega >= 1.7 || omega <= 1.9)), "wrong input for omega: " << omega);
    CHECK_MSG((norm > 0), "wrong input for normalizationfrequency: " << norm);
    CHECK_MSG((checkf > 0), "wrong input for checkfrequency: " << checkf);

    PROG("construct solver manually");
    // set solver values
    imax = i_max;
    jmax = j_max;
    smax = itermax;

    eps = ep;
    om = omega;
    normfreq = norm;
    checkfreq = checkf;
}

// Constructor to create a SORSolver from a parsed configuration file
SORSolver::SORSolver(const FileReader &configuration)
    : normfreq(1), checkfreq(1)
{
    PROG("construct solver with a file");
    // set solver values
    imax = configuration.getIntParameter("imax");
    CHECK_MSG((imax >= 0), "wrong input for imax: " << imax);
    jmax = configuration.getIntParameter("jmax");
    CHECK_MSG((jmax >= 0), "wrong input for jmax: " << jmax);
    smax = configuration.getIntParameter("itermax");
    CHECK_MSG((smax >= 0), "wrong input for itermax: " << smax);

    eps = configuration.getRealParameter("eps");
    CHECK_MSG((eps >= 0), "wrong input for eps: " << eps);
    om = configuration.getRealParameter("omg");
    CHECK_MSG((om == 1 || (om >= 1.7 || om <= 1.9)), "wrong input for omega: " << om);

    normfreq = configuration.getIntParameter("normalizationfrequency");
    if (normfreq == 0)
        normfreq = 1;
    CHECK_MSG((normfreq > 0), "wrong input for normalizationfrequency: " << normfreq);
    checkfreq = configuration.getIntParameter("checkfrequency");
    if (checkfreq == 0)
        checkfreq = 1;
    CHECK_MSG((checkfreq > 0), "wrong input for checkfrequency: " << checkfreq);

}


//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================

// solve the pressure equation on the staggered grid
bool SORSolver::solve(StaggeredGrid &grid)
{
    std::cout << "\n";
    PROG("check existence of solution");
    // check existence of solution
    real sum = 0;

    for (int i = 0; i < imax; ++i)
    {
        for (int j = 0; j < jmax; ++j)
        {

            if (grid.isFluid(i + 1, j + 1))
                sum += grid.rhs()(i, j);
        }
    }

    if (fabs(sum) > 0.000000001)
    {
        WARN("A solution doesn't exist. Instability!");

        real rhsMean = sum / grid.getNumFluid();
        // substract the mean
        for (int i = 0; i < imax; ++i)
        {
            for (int j = 0; j < jmax; ++j)
            {

                if (grid.isFluid(i + 1, j + 1))
                    grid.rhs()(i, j) -= rhsMean;
            }
        }
    }
    else
    {
        PROG("a solution exists");
    }

    // declaration of actual steps and residue
    int step = 0;
    real residue = 0.0;

    // calculate and check initial residue
    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
        {

            if (grid.isFluid(i, j))
                residue += (grid.rhs()(i - 1, j - 1) - grid.p()(i, j)) * (grid.rhs()(i - 1, j - 1) - grid.p()(i, j));
        }
    }
    residue = sqrt(residue / grid.getNumFluid());
    if (residue <= eps)
    {
        std::cout << "total iteration steps = " << step << "\n" << std::endl;

        return true; // stop if initial residue is ok
    }

    PROG("start iteration");

    real oneDXsq = 1 / (grid.dx() * grid.dx());
    real oneDYsq = 1 / (grid.dy() * grid.dy());
    real oneTerm = 1 / (2 * oneDXsq + 2 * oneDYsq);

    // iterate
    while (step <= smax)      // check actual steps
    {

        residue = 0.0;
        for (int i = 1; i <= imax; ++i)
        {
            for (int j = 1; j <= jmax; ++j)
            {

                if (grid.isFluid(i, j))
                {
                    // set new border values
                    if (i == 1)   // left
                        grid.p()(i - 1, j) = grid.p()(i, j);

                    if (i == imax)   // right
                        grid.p()(i + 1, j) = grid.p()(i, j);

                    if (j == 1)   // down
                        grid.p()(i, j - 1) = grid.p()(i, j);

                    if (j == jmax)   // up
                        grid.p()(i, j + 1) = grid.p()(i, j);

                    int number_of_empty_neighbour = 0 ;

                    // Count the number of empty neighbour cell
                    if (grid.isEmpty(i + 1, j))
                        ++ number_of_empty_neighbour ;
                    if (grid.isEmpty(i - 1, j))
                        ++ number_of_empty_neighbour ;
                    if (grid.isEmpty(i, j + 1))
                        ++ number_of_empty_neighbour ;
                    if (grid.isEmpty(i, j - 1))
                        ++ number_of_empty_neighbour ;

                    if (number_of_empty_neighbour == 0)
                    {
                        // calculate p
                        real oP = grid.p()(i, j); // old p (for residue)

                        real old = (1 - om) * grid.p()(i, j);
                        real factor = om * oneTerm;
                        real upgr1 = (grid.p(i + 1, j, WEST) + grid.p(i - 1, j, EAST)) * oneDXsq;
                        real upgr2 = (grid.p(i, j + 1, SOUTH) + grid.p(i, j - 1, NORTH)) * oneDYsq;
                        grid.p()(i, j) = old + factor * (upgr1 + upgr2 - grid.rhs()(i - 1, j - 1));   // actual iteration step

                        // calculate new "residue"
                        residue += (grid.p()(i, j) - oP) * (grid.p()(i, j) - oP);
                    }
                }
            }
        }
        // check new residue after checkfreq steps
        if (step % checkfreq == 0)
        {
            if (residue <= eps)
            {
                PROG("iteration was successful");
                std::cout << "total iteration steps = " << step << "\n" << std::endl;

                return true; // stop if new residue is ok
            }
        }

        // one step done
        step++;
        std::cout << "current iteration step is " << step << ": current residue is " << residue << std::endl;

    }

    PROG("iteration failed");
    std::cout << "total iteration steps = " << step << "\n" << std::endl;

    return false;
}
