#ifndef VTKFILEWRITER_HH
#define VTKFILEWRITER_HH


#include "StaggeredGrid.hh"
#include "ParticleTracer.hh"



//*******************************************************************************************************************
/*! Writes StaggeredGrid as vtk file

- vtk files can for example be opened with Paraview ( http://www.paraview.org/ )
- writes out pressure and/or velocity

- Usage:
\code
VTKWriter vtkWriter ( myGrid, "lidDrivenCavity", true, true );
// for each timestep:
vtkWriter.write();
\endcode
This creates on file per timestep: "lidDrivenCavity_0001.vtk", ""lidDrivenCavity_0002.vtk" ...

*/
//*******************************************************************************************************************
class VTKWriter
{

public:

    VTKWriter(const std::string &basename);

    void write(const StaggeredGrid &grid, const ParticleTracer *tracer);

private:
    std::string baseName_;

    int counter_;
};

#endif
