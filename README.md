# RMT3D: Reference Map Technique 3D
Authors: Yuexia Luna Lin, Nicholoas J. Derr, Chris H. Rycroft
RMT3D simulates 3D incompressible fluid-structure interaction with solids that undergoes large deformation.
An incompressible neo-Hookean material model is used.

## Software/library requirements:

- GNUMake
- GNUGCC compiler
- <a href="https://github.com/chr1shr/pgmg" >PGMG library </a>
- Perl         (For image processing only)
- Povray 3.7   (For image processing only)

## Get started:
1. Download or <code> git clone </code> RMT3D repository.
2. Download or <code> git clone </code> PGMG repository. For convenience, put PGMG/ and RMT3D/ on the same file system level.
3. Create a config.mk file to configure which compiler to use and basic compiler and linker flags. Sample config.mk files are provided for Mac and Linux in the make_configs/ folder.
4. Change into Execs/ folder and open Make.defs. Modify the line <code> PGMG_DIR = </code> to point to the path to PGMG/ repo. For example, <pre> PGMG_DIR = /Users/gitcodes/PGMG/ </pre> Notice that the library name is _libmg3d_.
6. After these changes, typing <code> make </code> in the commandline in Execs/ directory will build an MPI-enable application, with executable run_sim.
7. To run the application, a config file (must have file extension .cfg) for the simulation must be provided as command line argument. To run it in serially,  <pre> ./run_sim mySim.cfg </pre> To run it with N processes, <pre> mpirun -np N run_sim mySim.cfg </pre>
8. Sample simulation config files can be found in sim_configs/ directory. Examples include
 * a full fluid test case with manufactured solution
 * a full solid test case with solid shear wave
 * a lid-driven cavity flow (full fluid)
 * a lid-driven cavity flow with a sphere
 * a pre-stretched sphere relaxing in fluid
9. If POV-Ray is installed and the necessary output files (contours and tracers) are available, the perl script pov-movie.pl can be used to render the output files into 3D snapshots. For detailed usage, see the output of <pre> perl pov-movie.pl -h </pre>

More complete documentation is under development.
In the meantime, for any questions, feel free to contact the authors of this repository.

## Acknowledgement
This work has been partially supported by the Applied Mathematics Program of the U.S. DOE Office of Science Advanced Scientific Computing Research under contract number DE-AC02-05CH11231, the Department of Energy Computational Science Graduate Fellowship, the Department of Defense NDSEG Fellowship, and the Harvard NSF-Simons Center Quantitative Biology Initiative student fellowship, supported by the National Science Foundation grant DMS-1764269.

## References
Yuexia Luna Lin, Nicholas J. Derr, and Chris H. Rycroft, _Eulerian simulation of complex suspensions and biolocomotion in three dimensions_, <a href="https://arxiv.org/abs/2104.00095">arXiv:2104.00095 [physics.flu-dyn]</a>.

Chris H. Rycroft _et al._, _Reference map technique for incompressible fluid-structure interaction_, Journal of Fluid Mechanics, 898, A9 (2020). <a href= "https://doi.org/10.1017/jfm.2020.353"> doi:10.1017/jfm.2020.353 </a>


