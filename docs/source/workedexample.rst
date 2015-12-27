Parameterization
================

An example is placed in *example/jmak/SCORE.Input.TemplateIron.uds*.

Simulation
==========

The sequential (MPI) simulation is issued by the following command line call::

   mpiexec -n 1 score SCORE.Input.TemplateIron.uds 1

Output
======

The *example/jmak* folder provides a simple recrystallization case study of 1000 nuclei growing into a deformed monocrystal.
The output contains several files, all of which have a prefix **SCORE.<SimulationID>** to make distinguishing results simple.
Most output is plain ASCII formatted such that it can be directly imported in post-processing tools.

Rendering
=========

 .. figure:: ../images/ImportParaview_05.jpg
   :scale: 50%
   :align: center
