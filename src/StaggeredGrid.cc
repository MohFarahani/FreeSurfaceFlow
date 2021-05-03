
#include "StaggeredGrid.hh"

StaggeredGrid::StaggeredGrid(int xxSize, int yySize, real ddx, real ddy)
{
    CHECK_MSG((xxSize >= 0), "wrong input for xSize: " << xxSize);
    CHECK_MSG((yySize >= 0), "wrong input for ySize: " << yySize);
    CHECK_MSG((ddx >= 0), "wrong input for dx: " << ddx);
    CHECK_MSG((ddy >= 0), "wrong input for dy: " << ddy);
    PROG("construct grid manually");
    // initilize grid values
    dx_ = ddx;
    dy_ = ddy;
    xSize_ = xxSize;
    ySize_ = yySize;
    imax_ = (int)(xSize_ / dx_);
    jmax_ = (int)(ySize_ / dy_);
    ppc_ = 9;

    Array pp(imax_ + 2, jmax_ + 2);
    Array rhss(imax_, jmax_);
    Array uu(imax_ + 1, jmax_ + 2);
    Array vv(imax_ + 2, jmax_ + 1);
    Array ff(imax_ + 1, jmax_);
    Array gg(imax_, jmax_ + 1);
    FlagArray ob(imax_ + 2, jmax_ + 2);

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    fluid cell (center):    ob(i,j) =  2
    //    empty cell (center):    ob(i,j) =  4

    ob.fill(FLUID); 
    obs_ = ob;

}

StaggeredGrid::StaggeredGrid(const FileReader &configuration)
{
    PROG("construct grid with file");
    // initilize grid values

    real xl = configuration.getRealParameter("xlength");
    CHECK_MSG((xl >= 0), "wrong input for xlength: " << xl);
    real yl = configuration.getRealParameter("ylength");
    CHECK_MSG((yl >= 0), "wrong input for ylength: " << yl);
    int iimax = configuration.getIntParameter("imax");
    CHECK_MSG((iimax >= 0), "wrong input for imax: " << iimax);
    int jjmax = configuration.getIntParameter("jmax");
    CHECK_MSG((jjmax >= 0), "wrong input for jmax: " << jjmax);
    int pppc = configuration.getIntParameter("ppc");
    if (pppc <= 0)
        pppc = 9;

    dx_ = xl / iimax;
    dy_ = yl / jjmax;
    xSize_ = xl;
    ySize_ = yl;
    imax_ = iimax;
    jmax_ = jjmax;
    ppc_ = pppc;

    Array pp(imax_ + 2, jmax_ + 2);
    Array rhss(imax_, jmax_);
    Array uu(imax_ + 1, jmax_ + 2);
    Array vv(imax_ + 2, jmax_ + 1);
    Array ff(imax_ + 1, jmax_);
    Array gg(imax_, jmax_ + 1);
    FlagArray ob(imax_ + 2, jmax_ + 2);

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    fluid cell (center):    ob(i,j) =  2
    //    empty cell (center):    ob(i,j) =  4

    ob.fill(FLUID);
    obs_ = ob;

}

void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
    // define corners in grid
    int startX = std::min((int)(x1 / dx_), (int)(x2 / dx_));
    int endX = std::max((int)(x1 / dx_), (int)(x2 / dx_));
    int startY = std::min((int)(y1 / dy_), (int)(y2 / dy_));
    int endY = std::max((int)(y1 / dy_), (int)(y2 / dy_));

    for (int i = startX; i <= endX; ++i)
    {
        for (int j = startY; j <= endY; ++j)
        {
            //PROG("in for");
            setCellToObstacle(i, j); // set cell to obstacle
        }
    }

}

void StaggeredGrid::createCircle(real x, real y, real r)
{
    real point = 0.0;
    for (int i = 0; i < obs_.getSize(0); ++i)
    {
        for (int j = 0; j < obs_.getSize(1); ++j)
        {

            // compute vector to the actual point
            point = (x - i * dx_) * (x - i * dx_) + (y - j * dy_) * (y - j * dy_);

            if (point <= r * r) // check if point is in the circle
                setCellToObstacle(i, j); // set cell to obstacle
        }
    }

}

void StaggeredGrid::createPng(const std::string &pngFilename)
{
    // create png file
    int n = obs_.getSize(0);
    int m = obs_.getSize(1);
    GrayScaleImage image(pngFilename, n, m);

    // set size
    image = image.getResizedImage(n, m);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (isFluid(i, j))    // set cell to fluid
            {
                image(i, j) = std::numeric_limits<unsigned char>::max();
            }
            else     // set cell to obstacle
            {
                image(i, j) = 0;
            }
        }
    }

    // save image
    image.save(pngFilename);

}

void StaggeredGrid::readPng(const std::string &pngFilename)
{
    // load png file
    GrayScaleImage image(pngFilename);

    // get size
    int n = image.width();
    int m = image.height();

    // create flag array
    FlagArray ob(n, m);
    ob.fill(FLUID);

    // assign to obstacle array
    obs_ = ob;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (image.getElement(i, j) == 0)    // set cell to obstacle
            {
                setCellToObstacle(i, j);
            }
        }
    }

}

void StaggeredGrid::setCellToFluid(int x, int y)
{
    obs_(x, y) = FLUID;
}

void StaggeredGrid::setCellToEmpty(int x, int y)
{
    obs_(x, y) = EMPTY;
}

void StaggeredGrid::refreshEmpty()
{
    for (int i = 1; i <= imax_; ++i)
    {
        for (int j = 1; j <= jmax_; ++j)
        {
            if(!isEmpty(i,j)) continue;

            p_(i,j) = 0;

            if(!isFluid(i+1,j)) {
                u_(i, j) = 0;
            }

            if(!isFluid(i,j+1)) {
                v_(i, j) = 0;
            }
        }
    }
}

void StaggeredGrid::setCellToObstacle(int x, int y)
{
    obs_(x, y) = OBS; // set cell to obstacle

    // it is okay to set the velocities of obstacles right here
    // because they are never updated or moved
    // this IS NOT the case for empty cells which change constantly
    if (x < u_.getSize(0) && y < u_.getSize(1))
        u_(x, y) = 0;
    if (x < v_.getSize(0) && y < v_.getSize(1))
        v_(x, y) = 0;

    //     // configurate neighbors
    //     if ( (x - 1) >= 0 ) // east
    //         obs_(x - 1, y) = obs_(x - 1, y) | OBSWEST;
    //     if ( (x + 1) < obs_.getSize(0) ) // west
    //         obs_(x + 1, y) = obs_(x + 1, y) | OBSEAST;
    //     if ( (y - 1) >= 0 ) // south
    //         obs_(x, y - 1) = obs_(x, y - 1) | OBSNORTH;
    //     if ( (y + 1) < obs_.getSize(1) ) // north
    //         obs_(x, y + 1) = obs_(x, y + 1) | OBSSOUTH;
}
