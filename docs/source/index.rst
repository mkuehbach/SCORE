.. Documentation to the SCORE
.. ==========================

| SCORE is an MPI-parallelized 3D cellular automaton for the simulation of microstructure evolution during primary recrystallization in polycrystals. 
| The source code is developed by Markus Kühbach at the Institute of Physical Metallurgy and Metal Physics with additional contributions from Luis Antonio Barrales-Mora.

 .. figure:: ../images/SCORELogo.jpg
   :scale: 50%
   :align: left
   :target: http://www.imm.rwth-aachen.de/index.php?id=88&L=1
   
|
|
|
|
|
|
|
|
|

1. Getting started
^^^^^^^^^^^^^^^^^^

.. .. figure:: ../images/VoxelizedGrainsGeneric.jpg
   :scale: 40%
   :align: center

| The compilation of SCORE utilizes standard Linux tools.

.. toctree::
   :maxdepth: 2
     
   setup
   scorephysics
   
   
2. Parameterization
^^^^^^^^^^^^^^^^^^^

Parameterizing and setting up advanced recrystallization models can appear cumbersome at first glance owing to the variety of settings and input parameter. 
In order to resolve many of the questions, please find in the following a detailed description of all input parameter to the SCORE model.
The structure of the **.uds** file identifies where data should be read as strings ( **s** ), integers ( **i** ) or double precision floating point numbers ( **f** ).
In general, the parameter values are categorized loosely into these groups.

.. toctree::
   :maxdepth: 2
   
   parmconfig
   
   
   
3. Simulations
^^^^^^^^^^^^^^
 
A simulation is issued by the following command line call::

   mpiexec -n <proc> <scorempi> <Input.uds> <SimulationID>
   
*<scorempi>* denotes the name of the binary file. In its current version two input arguments are necessary:
 | <Input.uds> 
 |    this is the settings file of *uds* extension, relative path references are utilized by the program
 | <SimulationID> 
 |    this is a positive integer with which all simulation result can be uniquely labeled
 | **Mind that executing multiple simulations with the same ID in the same directory overwrites unpredictively previously obtained results without prompting!**

 

4. Model output
^^^^^^^^^^^^^^^

 .. figure:: ../images/GenericResultSingleKinetics.jpg
   :scale: 40%
   :align: center
   

Kinetics, macrotexture and grain size distribution are the key state parameters of interest. 

.. toctree::
   :maxdepth: 2
   
   theoutput
   
   
5. Visualize
^^^^^^^^^^^^

 .. figure:: ../images/TWIP3DCAExample.jpg
   :scale: 50%
   :align: center

| The 3D voxelized microstructure snapshots that become generated read as plain implicit 3D arrays. 
| As such, they can readily be imported into visualization software, the details how to do so are described below:

.. toctree::
   :maxdepth: 2

   score2pview

   
.. 6. Worked-out example
.. ^^^^^^^^^^^^^^^^^^^^^
.. 
.. In a nutshell, all this is exemplified on a small JMAK-type primary recrystallization simulation that can be found in .. the *example/* folder (SCORE.Input.1.uds).

.. .. toctree::
..    :maxdepth: 2
..    
..    workedexample
   
   
References
^^^^^^^^^^

**This model:**
 | Kühbach, M., Gottstein, G., Barrales-Mora, L.A.
 | A Statistical Ensemble Cellular Automaton Microstructure Model for Primary Recrystallization
 | submitted to "Acta Materialia" on Dec, 27, 2015
.. http://dx.doi.org/doi:10.1016/j.actamat.2011.07.003 
 | http://github.com/mkuehbach/SCORE

 
.. toctree::
  :maxdepth: 2
  
  references
  
 >
**Funding:**
 | The authors gratefully acknowledge the support from the DFG in the frame of the Reinhart Koselleck project (GO 335/44-1).
 | Furthermore, we acknowledge the support from the FZJuelich and RWTH Aachen University within the project JARAHPC projects. 
 
 
Version history
^^^^^^^^^^^^^^^

**v1.0** first version, 
 |     with deformed grain boundary detection and nucleation model, Poisson-Voronoi structure, Poissonian nucleation
 |     site-saturated and time-dependent nucleation, MPI-parallelized, summary statistics and basic voxel structures


MPI/OpenMP extension
^^^^^^^^^^^^^^^^^^^^

In addition the code has been modified for thread-parallelism with OpenMP, and the source code is made available in the **SCORE_MPIOMP_20151224.zip** file that resides in the **hybrid** folder. The user is motivated to extract the archive in a separate folder and to build the program as it was described already for the MPI-only version. 

This hybrid SCORE model works as the MPI-parallelized version, and enabled a reduction of the execution time by a factor of 9, at most, as it was benchmarked on a Bullx X5675 Xeon two socket cluster blade system with exclusive usage of the 12 physical cores (OMP_NUM_THREADS=12) and explicitly bound the threads to the cores (KMP_AFFINITY). It is possible already with this model to handle Poissonian nucleation in cuboidal deformed grain structures.

However, **not all internal functions have already been thread-parallelized, such as the grain boundary identification and grain boundary nucleation, and hence these functionalities are claimed as disfunctional at the moment**.

.. toctree::
  :maxdepth: 2
  
  openmp

 
   
Licence
^^^^^^^

The project is licenced under the GNU v2.0


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

Questions, contributions
^^^^^^^^^^^^^^^^^^^^^^^^

Just let me know or contact *markus.kuehbach@rwth-aachen.de*
