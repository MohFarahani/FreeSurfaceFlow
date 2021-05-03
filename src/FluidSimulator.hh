#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "ParticleTracer.hh"

#include "VTKWriter.hh"

#include <cmath>


class FluidSimulator
{
public:
    FluidSimulator(const FileReader &conf);

    /// Simulates a given time-length
    void simulate(real duration);
    void simulateTimeStepCount(unsigned int nrOfTimeSteps);


    // Getter functions for the internally stored StaggeredGrid
    StaggeredGrid &grid()
    {
        return grid_;
    };
    const StaggeredGrid &grid() const
    {
        return grid_;
    };

    // Getter functions for the internally stored SORSolver
    SORSolver &solv()
    {
        return solver_;
    };
    const SORSolver &solv() const
    {
        return solver_;
    };

    ParticleTracer &tracer()
    {
        return particle_tracer_;
    };
    const ParticleTracer &tracer() const
    {
        return particle_tracer_;
    };


    void test();
    void normalization();
    void dcavity_test();
    void testFG();

    // only public for tests
    void set_UVP_surface(const real &dt, bool compP);

private:
    // helper functions (derivatives)
    real dxuu(int i, int j), dyuv(int i, int j), ddxu(int i, int j), ddyu(int i, int j);
    real dyvv(int i, int j), dxuv(int i, int j), ddxv(int i, int j), ddyv(int i, int j);

    void computeFG();
    void composeRHS();
    void updateVelocities();
    void determineNextDT(real const &limit);
    void refreshBoundaries();

    // compute surface boundary
    void set_UVP_surface(int i, int j , const real &dt, bool compP);
    void one_empty_neighbour(int i , int j , const real &dt, bool compP) ;
    void two_empty_neighbour(int i , int j , const real &dt, bool compP) ;
    void three_empty_neighbour(int i , int j , const real &dt, bool compP) ;
    void four_empty_neighbour(int i , int j , const real &dt, bool compP) ;

    // needed values
    real safetyfac_, gamma_, Re_, gx_, gy_, dt_, vel_N, vel_S, vel_E, vel_W, uInit_, vInit_, pInit_;
    real rectX_, rectXX_, rectY_, rectYY_, circX_, circY_, circR_;
    BCTYPE cond_N, cond_S, cond_E, cond_W;
    unsigned int timeStepNr, normfreq, outPutInt;
    std::string name_;
    int imax, jmax;
    real rectX1_particle_ , rectX2_particle_ , rectY1_particle_ , rectY2_particle_, circX_particle_, circY_particle_, circR_particle_;

protected:
    StaggeredGrid grid_;   //< grid
    SORSolver solver_;     //< solver
    ParticleTracer particle_tracer_; //< particle tracer

};



#endif
