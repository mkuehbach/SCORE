| It is practical to distinguish the output in two categories.
| The first considers the simulated state variable values and a depiction of the microstructure. 
| The second considers information of the model performance.

All output from a simulation is written in the output folder using relative addressing.
The *example/jmak* folder provides a simple recrystallization case study of 1000 nuclei growing into a deformed monocrystal.
The output contains several files, all of which have a prefix **SCORE.<SimulationID>** to make distinguishing results simple.
Most output is plain ASCII formatted in such a way that it can be directly imported in post-processing, and spreadsheet tools.

Recrystallization kinetics
==========================

The **Rediscretized.Kinetics.csv** file contains a table with the typical time, recrystallized fraction as well as the Avrami plot logarithmic linearization.
The kinetics of each automaton can be obtained by switching on the **OutputRXFrontStats** with the first column giving the individual consumption of each domain and the further columns provide details about the microstructural path an the memory utilization of the interface cell lists.


Macrotexture evolution
======================

The **Rediscretized.Macrotexture.csv** file contains a table with the evolution of the macrotexture components as they were defined as **IdealOrientations**.
Counts are listed in cells summarized and rediscretized over all domains with the total amount of cells in the ensemble to the right.

Grain size distribution
=======================

The **Rediscretized.FinalGrainSize.csv** file contains a table with the final size of all grains aggregated over all domains. 
The header line provides basic descriptive statistics of the whole ensemble.
Beneath, a table contains a sorted list of all grains with their size, the domain they grew, in their orientation class as an index of the
n-th **IdealOrientation**, the normalized cumulative grain size distribution to inspect tail departures, and the incubation time after which they were instantiated in the simulation. Noteworthy, this file contains only nuclei that participated in the transformation.

2D sections
===========

A section of the microstructure at various recrystallized volume fractions can be printed. For this logpoints have to be defined in the input file block named **RenderMicrostructure**. Exemplarily, if one desires to output the microstructure at 25% and 76.8% recrystallized volume fraction, two lines should be placed in this block.

| "X"		0.25
| "X"		0.768

Lode Vandevenne's open source lodePNG (https://github.com/lvandeve/lodepng/) tool is utilized to render indexed RGB **png** files that are stored in the output directory. The coordinate origin is the lower left corner. The coordinate system identifies x as the first coordinate, y as the second and z as the third. This is to of relevance for all implicit addressing in the source code. By default the whole domain is sectioned at the z coordinate equal to 0.5.

3D data
======= 

It is also possible to render three-dimensionally the structure. **Mind, though, that this equals a file size of 4 byte per cell, i.e. a simulation with only 1000 x 1000 x 1000 already requires 4 GB of disk space.**
However, it is possible to restrict which simulation domains are rendered by indenting explicitly in the input file
the respective automaton IDs under **RenderMicrostructure**.

These files are stored as binary **raw** files in the working directory. At the moment there is no header associated with the file so remembering their size is important.

**For both 2D and 3D data two coloring modes are currently implemented.**
The first is coloring according to the grain ID. The second is to colorize according to inverse polefigure orientation mapping in the stereographic standard triangle of under fcc crystal symmetry and triclinic sample symmetry.


Single grain growth charts
==========================

The **SingleGrainData** sheet compiles the detailed evolution, nucleation sites, orientation, state of the deformed grains, automaton-resolved to allow insights into the evolution of particular grains, their incubation time, and much more. Furthermore, but currently commented out is an option to export these data by plain MPI I/O for post-processing purposes.


Performance
===========

First the **ProfilingLog** csv file present a summary of all automata that were executed. Along with MPI_Wtime-ings in front and past barriers,
in conjunction with a detailed information about how many nuclei in which orientation into which deformation structure,
this file is the first address for a statistical interpretation of the data.

In particular for developers the **RXFrontsStats** csv file provides detailed information about the course of microstructure evolution as well as
microstructural path parameter, discrete interfacial area the status of the recrystallization front memory management containers, etc.

The developmental hybrid (MPI/OpenMP) parallelized extension of the model adds a further output file **OMPThreadProfiling** which details how much time the OpenMP threads spent in each section of the code and synchronization.

