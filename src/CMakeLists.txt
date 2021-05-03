
include_directories( . )

add_executable( nusif FileReader.cc Debug.cc Array.cc SORSolver.cc StaggeredGrid.cc FluidSimulator.cc Particle.cc ParticleTracer.cc VTKWriter.cc GrayScaleImage.cc lodepng.cc main.cc )
