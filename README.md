# SCORE
MPI-parallelized cellular automaton ensemble model for the simulation of primary recrystallization phenomena in 3D

A detailed documentation of the model and the program is available under:
http://score.readthedocs.org/en/master/

The authors of SCORE gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft (DFG) within the "Reinhart Koselleck-Project" (GO 335/44-1) and the granting of computing time which was provided by the RWTH Aachen University JARAHPC research alliance (JARA0076).

The authors of SCORE gratefully acknowledge the financial support from the Max-Planck-Gesellschaft through the provisioning of computing time grants
in the frame of the BiGmax, the Max-Planck-Society's Research Network on Big-Data-Driven Materials Science.


Latest bugfixes:
- 2016/01/13: uncomment the line #define DEBUG in src/SCORE_Io.h when utilizing the GNU compiler
- 2019/02/23: we implemented a right-handed coordinate system and made modifications to use the code
for the following paper. M. Diehl and M. Kuehbach, Coupled experimental-computational analysis of primary static recrystallization in low carbon steel",
Modelling and Simulation in Materials Science and Engineering, 28, 2020, http://doi.org/10.1088/1361-651X/ab51bd
