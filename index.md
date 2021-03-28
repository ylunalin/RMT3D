# RMT3D
Authors: Yuexia Luna Lin, Nicholoas J. Derr, Chris H. Rycroft

The code simulates 3D incompressible fluid-structure interaction with neo-Hookean materials via the Reference Map Technique.

Software/library requirements:

- GNUMake
- GNUGCC compiler
- MG3D library (link)
- Perl         (For image processing only)
- Povray 3.7   (For image processing only)

Get started:
1. Download or <code> git clone </code> RMT3D repository.
2. Download or <code> git clone </code> MG3D repository. For convenience, put MG3D/ and RMT3D/ on the same file system level.
3. Create a config.mk file to configure which compiler to use and basic compiler and linker flags. Sample config.mk files are provided for Mac and Linux in the make_configs/ folder.
4. Change into Execs/ folder and open Make.defs. Modify the line <code> MG3D_DIR = </code> to point to the path to MG3D/ repo. For example, <pre> MG3D_DIR = /Users/gitcodes/MG3D/ </pre>
5. After these changes, typing <code> make </code> in the commandline in Execs/ directory will build an MPI-enable application, with executable run_sim.
6. To run the application, a config file (must have file extension .cfg) for the simulation must be provided as command line argument. To run it in serially,  <pre> ./run_sim mySim.cfg </pre> To run it with N processes, <pre> mpirun -np N run_sim mySim.cfg </pre>
7. Sample simulation config files can be found in sim_configs/ directory. Examples include
 * a full fluid test case with manufactured solution
 * a full solid test case with solid shear wave
 * a lid-driven cavity flow (full fluid)
 * a lid-driven cavity flow with a sphere
 * a pre-stretched sphere relaxing in fluid
8. If POV-Ray is installed and the necessary output files (contours and tracers) are available, the perl script pov-movie.pl can be used to render the output files into 3D snapshots. For detailed usage, see the output of <pre> perl pov-movie.pl -h </pre>

More complete documentation is under development.
In the meantime, for any questions, feel free to contact the authors of this repository.
