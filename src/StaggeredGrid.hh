#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH

#include <vector>

#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"
#include "GrayScaleImage.hh"
#include "Particle.hh"


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:

    // Constructors to manually create staggered grid
    StaggeredGrid(int xSize, int ySize, real ddx, real ddy);

    // Constructor to create a staggered grid from a parsed configuration file
    StaggeredGrid(const FileReader &configuration);

    // Assignment Operator
    inline StaggeredGrid &operator = (const StaggeredGrid &s);

    // Getters / Setters for member variables
    Array &p()
    {
        return p_;
    }
    Array &rhs()
    {
        return rhs_;
    }
    Array &u()
    {
        return u_;
    }
    Array &v()
    {
        return v_;
    }
    Array &f()
    {
        return f_;
    }
    Array &g()
    {
        return g_;
    }
    FlagArray &obs()
    {
        return obs_;
    }

    const Array &p()   const
    {
        return p_;
    }
    const Array &rhs() const
    {
        return rhs_;
    }
    const Array &u()   const
    {
        return u_;
    }
    const Array &v()   const
    {
        return v_;
    }
    const Array &f()   const
    {
        return f_;
    }
    const Array &g()   const
    {
        return g_;
    }
    const FlagArray &obs() const
    {
        return obs_;
    }
    const std::vector<Particle> &particles() const
    {
        return particles_;
    }

    real dx() const
    {
        return dx_;
    }
    real dy() const
    {
        return dy_;
    }

    real xSize() const
    {
        return xSize_ / dx_;
    }
    real ySize() const
    {
        return ySize_ / dy_;
    }

    real xDomain() const
    {
        return xSize_;
    }

    real yDomain() const
    {
        return ySize_;
    }

    int imax() const
    {
        return imax_;
    }
    int jmax() const
    {
        return jmax_;
    }

    int ppc() const
    {
        if (ppc_ != 0)
            return ppc_;

        return 9;
    }

    inline bool isFluid(const int x, const int y);
    inline bool isEmpty(const int x, const int y);
    inline bool isObstacle(const int x, const int y);
    inline int getNumFluid();

    inline real u(const int x, const int y, Direction dir);
    inline real v(const int x, const int y, Direction dir);
    inline real p(const int x, const int y, Direction dir);

    void setCellToFluid(int x, int y);
    void setCellToEmpty(int x, int y);
    void refreshEmpty();
    void setCellToObstacle(int x, int y);
    void createRectangle(real x1, real y1, real x2, real y2);
    void createCircle(real x, real y, real r);
    void createPng(const std::string &pngFilename);
    void readPng(const std::string &pngFilename);
    void markCells();

    //    // Interpolated function ( bilinear interpolation )
    //    real u_inter ( real x , real y ) ;
    //    real v_inter ( real x , real y ) ;
protected:
    Array p_;   //< pressure field
    Array rhs_; //< right hand side of the pressure equation
    Array u_;   //< velocity in x direction
    Array v_;   //< velocity in y direction
    Array f_;   //< f
    Array g_;   //< g
    FlagArray obs_; //< obstacle array with flags (1 for fluid cell, 0 for obstacle cell)

    real dx_;   //< distance between two grid points in x direction
    real dy_;   //< distance between two grid points in y direction
    real xSize_;
    real ySize_;
    int imax_, jmax_;
    int ppc_;

    std::vector<Particle> particles_;
};



inline bool StaggeredGrid::isFluid(const int x, const int y)
{
    return (obs_(x, y) & FLUID);
}

inline bool StaggeredGrid::isEmpty(const int x, const int y)
{
    return (obs_(x, y) & EMPTY);
}

inline bool StaggeredGrid::isObstacle(const int x, const int y)
{
    return (obs_(x, y) & OBS);
}


inline int StaggeredGrid::getNumFluid()
{
    int sum = 0;
    for (int i = 1; i <= imax_; ++i)
    {
        for (int j = 1; j <= jmax_; ++j)
        {
            int number_of_empty_neighbour = 0 ;

            if (isFluid(i, j))
            {
                // Count the number of empty neighbour cell
                if (isEmpty(i + 1, j))
                    ++ number_of_empty_neighbour ;
                if (isEmpty(i - 1, j))
                    ++ number_of_empty_neighbour ;
                if (isEmpty(i, j + 1))
                    ++ number_of_empty_neighbour ;
                if (isEmpty(i, j - 1))
                    ++ number_of_empty_neighbour ;

                if (number_of_empty_neighbour == 0)
                {
                    ++sum;
                }
            }
        }
    }

    return sum;
}

// inline real StaggeredGrid::u(const int x, const int y, Direction dir)
// {
//     if (isFluid(x, y) && isFluid(x + 1, y))
//         return u_(x, y);
//
//     if (dir == NORTH)
//     {
//         if (!isFluid(x, y) && !isFluid(x + 1, y))
//             return -u_(x, y + 1);
//     }
//
//     if (dir == SOUTH)
//     {
//         if (!isFluid(x, y) && !isFluid(x + 1, y))
//             return -u_(x, y - 1);
//     }
//
//     return 0;
// }
//
//
// inline real StaggeredGrid::v(const int x, const int y, Direction dir)
// {
//     if (isFluid(x, y) && isFluid(x, y + 1))
//         return v_(x, y);
//
//     if (dir == WEST)
//     {
//         if (!isFluid(x, y) && !isFluid(x, y + 1))
//             return -v_(x - 1, y);
//     }
//
//     if (dir == EAST)
//     {
//         if (!isFluid(x, y) && !isFluid(x, y + 1))
//             return -v_(x + 1, y);
//     }
//
//     return 0;
//
// }
//
//
// inline real StaggeredGrid::p(const int x, const int y, Direction dir)
// {
//     if (!isFluid(x, y))
//     {
//         if (dir == NORTH)
//         {
//             return p_(x, y + 1);
//
//         }
//         else if (dir == SOUTH)
//         {
//             return p_(x, y - 1);
//
//         }
//         else if (dir == WEST)
//         {
//             return p_(x - 1, y);
//
//         }
//         else
//         {
//             return p_(x + 1, y);
//         }
//     }
//
//     return p_(x, y);
// }

inline real StaggeredGrid::u(const int x, const int y, Direction dir)
{
    //  if (!( (x + 1 < u_.getSize(0)) && (x >= 0) && (y >= 0) && (y < u_.getSize(1)) ))
    //      return 0; // desert index out of bounds

    if (!isObstacle(x, y) && !isObstacle(x + 1, y))
        return u_(x, y);

    if (dir == NORTH /*&& (y + 1 < u_.getSize(1))*/) // desert index out of bounds
    {
        if (isObstacle(x, y) && isObstacle(x + 1, y))
            return -u_(x, y + 1);
    }

    if (dir == SOUTH /*&& (y - 1 >= 0)*/) // desert index out of bounds
    {
        if (isObstacle(x, y) && isObstacle(x + 1, y))
            return -u_(x, y - 1);
    }

    if (dir == DIAG)
    {
        //      if (y + 1 < u_.getSize(1)) // desert index out of bounds
        //      {
        if (isObstacle(x, y) && isObstacle(x + 1, y) && isFluid(x, y + 1))
            return -u_(x, y + 1);
        //      }

        //      if (y - 1 >= 0) // desert index out of bounds
        //      {
        if (isObstacle(x, y) && isObstacle(x + 1, y) && isFluid(x, y - 1))
            return -u_(x, y - 1);
        //      }
    }

    return 0;
}


inline real StaggeredGrid::v(const int x, const int y, Direction dir)
{
    //  if (!( (y + 1 < v_.getSize(1)) && (y >= 0) && (x >= 0) && (x < v_.getSize(0)) ))
    //      return 0; // desert index out of bounds

    if (!isObstacle(x, y) && !isObstacle(x, y + 1))
        return v_(x, y);

    if (dir == WEST /*&& (x - 1 >= 0)*/) // desert index out of bounds
    {
        if (isObstacle(x, y) && isObstacle(x, y + 1))
            return -v_(x - 1, y);
    }

    if (dir == EAST /*&& (x + 1 < v_.getSize(0))*/) // desert index out of bounds
    {
        if (isObstacle(x, y) && isObstacle(x, y + 1))
            return -v_(x + 1, y);
    }

    if (dir == DIAG)
    {
        //      if (x + 1 < v_.getSize(0)) // desert index out of bounds
        //      {
        if (isObstacle(x, y) && isObstacle(x, y + 1) && isFluid(x + 1, y))
            return -v_(x + 1, y);
        //      }
        //
        //      if (x - 1 >= 0) // desert index out of bounds
        //      {
        if (isObstacle(x, y) && isObstacle(x, y + 1) && isFluid(x - 1, y))
            return -v_(x - 1, y);
        //      }
    }

    return 0;

}


inline real StaggeredGrid::p(const int x, const int y, Direction dir)
{
    if (isObstacle(x, y))
    {
        if (dir == NORTH)
        {
            return p_(x, y + 1);

        }
        else if (dir == SOUTH)
        {
            return p_(x, y - 1);

        }
        else if (dir == WEST)
        {
            return p_(x - 1, y);

        }
        else if (dir == EAST)
        {
            return p_(x + 1, y);
        }
        else // diagonal
        {
            if (!isObstacle(x + 1, y) && !isObstacle(x, y + 1)) // NE boundary of obstacle
                return 0.5 * (p_(x + 1, y) + p_(x, y + 1));

            if (!isObstacle(x - 1, y) && !isObstacle(x, y + 1)) // NW boundary of obstacle
                return 0.5 * (p_(x - 1, y) + p_(x, y + 1));

            if (!isObstacle(x + 1, y) && !isObstacle(x, y - 1)) // SE boundary of obstacle
                return 0.5 * (p_(x + 1, y) + p_(x, y - 1));

            if (!isObstacle(x - 1, y) && !isObstacle(x, y + 1)) // SW boundary of obstacle
                return 0.5 * (p_(x - 1, y) + p_(x, y + 1));

            if (!isObstacle(x, y + 1)) // N boundary of obstacle
                return p_(x, y + 1);

            if (!isObstacle(x, y - 1)) // S boundary of obstacle
                return p_(x, y - 1);

            if (!isObstacle(x + 1, y)) // E boundary of obstacle
                return p_(x + 1, y);

            if (!isObstacle(x - 1, y)) // W boundary of obstacle
                return p_(x - 1, y);
        }
    }

    return p_(x, y);
}

#endif //STAGGERED_GRID_HH
