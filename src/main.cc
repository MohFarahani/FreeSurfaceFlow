
#include "SORSolver.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"
#include "FileReader.hh"



int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cerr << "No config file given" << std::endl;
        return EXIT_FAILURE;
    }

    FileReader confi;
    confi.registerStringParameter("name");

    //////////////////////////////////////////////////////// Initialization ///////////////////////////////////////////////////////
    confi.registerStringParameter("boundary_condition_N");
    confi.registerStringParameter("boundary_condition_S");
    confi.registerStringParameter("boundary_condition_E");
    confi.registerStringParameter("boundary_condition_W");

    confi.registerRealParameter("boundary_velocity_N");
    confi.registerRealParameter("boundary_velocity_S");
    confi.registerRealParameter("boundary_velocity_E");
    confi.registerRealParameter("boundary_velocity_W");

    confi.registerRealParameter("GX");
    confi.registerRealParameter("GY");

    confi.registerRealParameter("Re");

    confi.registerRealParameter("U_INIT");
    confi.registerRealParameter("V_INIT");
    confi.registerRealParameter("P_INIT");

    //////////////////////////////////////////////////////// Geometry Data ////////////////////////////////////////////////////////
    confi.registerRealParameter("xlength");
    confi.registerRealParameter("ylength");
    confi.registerIntParameter("imax");
    confi.registerIntParameter("jmax");
	
	confi.registerIntParameter("ppc");
	
    confi.registerRealParameter("RectangleParticleX1");
    confi.registerRealParameter("RectangleParticleX2");
    confi.registerRealParameter("RectangleParticleY1");
    confi.registerRealParameter("RectangleParticleY2");

    confi.registerRealParameter("CircleParticleX");
    confi.registerRealParameter("CircleParticleY");
    confi.registerRealParameter("CircleParticleR");

    confi.registerRealParameter("RectangleX1");
    confi.registerRealParameter("RectangleY1");
    confi.registerRealParameter("RectangleX2");
    confi.registerRealParameter("RectangleY2");

    confi.registerRealParameter("CircleX");
    confi.registerRealParameter("CircleY");
    confi.registerRealParameter("CircleR");

    //////////////////////////////////////////////////////// Time Data ////////////////////////////////////////////////////////////
    confi.registerRealParameter("dt");
    confi.registerIntParameter("timesteps");
    confi.registerRealParameter("safetyfactor");

    //////////////////////////////////////////////////////// Pressure Iteration Data //////////////////////////////////////////////
    confi.registerIntParameter("itermax");
    confi.registerRealParameter("eps");
    confi.registerRealParameter("omg");
    confi.registerRealParameter("gamma");
    confi.registerIntParameter("checkfrequency");
    confi.registerIntParameter("normalizationfrequency");

    //////////////////////////////////////////////////////// VTK Visualization Data ///////////////////////////////////////////////
    confi.registerIntParameter("outputinterval");



    //////////////////////////////////////////////////////// Test ////////////////////////////////////////////////////////////

    CHECK_MSG(confi.readFile(argv[1]), "Could not open file " << argv[1] << " which has to be in the current directory.");

    // create simulator
    FluidSimulator simulator(confi);

    // test
    simulator.test();

    std::cout << "Program end with: " << argv[1] << std::endl;

    return 0;
}

