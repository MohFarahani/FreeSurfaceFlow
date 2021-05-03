#include "VTKWriter.hh"
#include "Debug.hh"

#include <fstream>
#include <sstream>
#include <iomanip>


template<typename T> struct RealTypeToString         {};
template<>           struct RealTypeToString<float>
{
    static const char *str;
};
template<>           struct RealTypeToString<double>
{
    static const char *str;
};

const char *RealTypeToString<float>::str  = "float";
const char *RealTypeToString<double>::str = "double";



VTKWriter::VTKWriter(const std::string &basename)
    : baseName_(basename),
      counter_(0)
{
}

void VTKWriter::write(const StaggeredGrid &grid, const ParticleTracer *tracer)
{
    std::stringstream fileNameFluid;
    fileNameFluid << baseName_ << "_" <<  std::setw(4) << std::setfill('0') << counter_ << "_fluid.vtk";
    std::ofstream fileStreamFluid(fileNameFluid.str().c_str());

    fileStreamFluid << "# vtk DataFile Version 4.0\n";
    fileStreamFluid << "Nusif Fluid VTK output\n";
    fileStreamFluid << "ASCII\n";
    fileStreamFluid << "DATASET STRUCTURED_POINTS\n";

    fileStreamFluid << "DIMENSIONS " << grid.xSize() << " " << grid.ySize() << " 1\n";
    fileStreamFluid << "ORIGIN 0 0 0 \n";
    fileStreamFluid << "SPACING " << grid.dx() << " " << grid.dy() << " 1\n";
    fileStreamFluid << "POINT_DATA " << grid.xSize() * grid.ySize() << " \n" << std::endl;

    fileStreamFluid << "VECTORS velocity " << RealTypeToString<real>::str << "\n";

    for (int j = 0; j < grid.ySize(); ++j)
    {
        for (int i = 0; i < grid.xSize(); ++i)
        {
            const real u = 0.5 * (grid.u()(i + 1, j + 1)  + grid.u()(i, j + 1));
            const real v = 0.5 * (grid.v()(i + 1, j + 1)  + grid.v()(i + 1, j));
            fileStreamFluid << u << " " << v << " " << " 0\n";
        }
    }

    fileStreamFluid << "\n";

    fileStreamFluid << "SCALARS pressure " << RealTypeToString<real>::str << " 1\n";
    fileStreamFluid << "LOOKUP_TABLE default\n";

    for (int j = 0; j < grid.ySize(); ++j)
    {
        for (int i = 0; i < grid.xSize(); ++i)
        {
            fileStreamFluid << grid.p()(i + 1, j + 1) << "\n";
        }
    }

    fileStreamFluid << "SCALARS rhs " << RealTypeToString<real>::str << " 1\n";
    fileStreamFluid << "LOOKUP_TABLE default\n";

    for (int j = 0; j < grid.ySize(); ++j)
    {
        for (int i = 0; i < grid.xSize(); ++i)
        {
            fileStreamFluid << grid.rhs()(i, j) << "\n";
        }
    }

    if (tracer != NULL)
    {
        std::stringstream fileNameParticles;
        fileNameParticles << baseName_ << "_" <<  std::setw(4) << std::setfill('0') << counter_ << "_particles.vtk";
        std::ofstream fileStreamParticles(fileNameParticles.str().c_str());

        fileStreamParticles << "# vtk DataFile Version 4.0\n";
        fileStreamParticles << "Nusif Particles VTK output\n";
        fileStreamParticles << "ASCII\n";
        fileStreamParticles << "DATASET UNSTRUCTURED_GRID\n";

        fileStreamParticles << "POINTS " << tracer->particles().size() << " " << RealTypeToString<real>::str << std::endl;
        for (std::vector<Particle>::const_iterator p = tracer->particles().begin() ; p != tracer->particles().end(); ++p)
        {
            fileStreamParticles << std::fixed << std::setprecision(5) << p->x() - 0.5 * grid.dy() << " " << p->y() - 0.5 * grid.dy() << " 0\n";
        }
        fileStreamParticles << "CELLS 0 0" << std::endl;
        fileStreamParticles << "CELL_TYPES 0" << std::endl;
        fileStreamParticles << "POINT_DATA " << tracer->particles().size() << std::endl;
        fileStreamParticles << "SCALARS type int" << RealTypeToString<real>::str << std::endl;
        fileStreamParticles << "LOOKUP_TABLE default" << std::endl;
        for (std::vector<Particle>::const_iterator p = tracer->particles().begin() ; p != tracer->particles().end(); ++p)
        {
            fileStreamParticles << p->type() << std::endl;
        }

        fileStreamParticles << "\n";
    }

    ++counter_;
}

