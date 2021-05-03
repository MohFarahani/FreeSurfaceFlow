#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <iostream>
#include "Types.hh"

class Particle
{
public:

    Particle(real x, real y, int type);

    real &x()
    {
        return x_;
    }
    real &y()
    {
        return y_;
    }
    const real &x() const
    {
        return x_;
    }
    const real &y() const
    {
        return y_;
    }
    const int &type() const
    {
        return type_;
    }

    int getCellX(real dx);
    int getCellY(real dy);

    void setX(real x);
    void setY(real y);

private:
    real x_ ;   // "x" position of particle
    real y_ ;   // "y" position of particle

    int type_;  // encodes different particle types used for visualization
};

inline void Particle::setX(real xx)
{
    x_ = xx;
}
inline void Particle::setY(real yy)
{
    y_ = yy;
}

#endif //PARTICLE_HHH
