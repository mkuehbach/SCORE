Learn about the physical background to the model
================================================

 .. figure:: ../images/AA8079LCR706min300C.png
   :scale: 100%
   :align: center
   
Band contrast SEM/EBSD image of an electropolished AA8079 sample after 70% cold-rolling and isothermal annealing at 300 deg C for 6 minutes. 
.. The rolling direction is parallel to the scale bar. Blue denotes grains in cube orientation, high-angle grain boundaries drawn in black.

Recrystallization and recovery term a variety of deformation microstructure restoration processes during which the stored elastic energy in the deformed microstructure drives the elimination of excess dislocations and the migration of grain boundaries. Specifically the *static recrystallization* refers to a form of recrystallization which occurs after the deformation was applied. For microstructure evolution during the annealing of cold-rolled sheet material
this static recrystallization is of significant technical interest. The physical mechanisms and the strength of the driving forces distinguishes it from static recovery and other forms of grain growth.


**Scope of the model:**

The SCORE model is designed for the simulation of the growth of new nuclei under the effect of anisotropic
grain boundary mobility and heterogeneous deformation microstructures, i.e. for simulating *static recrystallization*.
Furthermore, SCORE encompasses approximate models for the nucleation stage. As with every model the nature of the processes is depicted in simplified manner by making the following assumptions:

* A microstructure volume comprises cubic voxel as the smallest discretization unit.
* The deformation microstructure can be idealized in contiguous regions with homogeneous properties (orientation, dislocation density).
* Therein nuclei appear which are strain- and excess dislocation free.

* Grain boundary migration follows classical rate theory according to Turnbull's model.
* Grain boundary migration is simplified as a flat boundary sweeping the deformation structure.
* The intrinsic grain boundary mobility is only dependent on the crystallographic disorientation of the adjoining grains.
* The deformation microstructure remains static in terms of orientation and boundaries content.
* However the (statistically stored) dislocation density is allowed to evolve. This is referred to as recovery modeling, for which approximate physical models exist.
* Recrystallizing nuclei do not rotate or change their orientation during growth.
* A potential growth inhibitating effect of second-phase particles is modeled in accordance with the Zener-Smith theory.


**Relevance of the model:**

These conditions apply well to many cases of cold-rolled aluminium, nickel, copper alloys and steels, but require that the driving force stemming for dislocations is at least an order of magnitude larger than that of curvature.
Under these conditions the SCORE model can predict the course of microstructure transformation in particular kinetics, indicate potential changes in the macrotexture, resolv spatially the impingement of individual nuclei and hence allows to extract very accurately the grain size distribution. The benefit of the model compared to similar approaches is
that it enables the simultaneous initialization and execution of multiple computational domains over which the results
can be averaged or analyzed in all details collectively. This renders the simulations statistical significant
and sufficiently discretized to minimize bias and to quantify uncertainty that persist to unknown extend in other recrystallization models.



**Further reading:**
 | Cotterill, P., Mould, P. R.
 | Recrystallization and Grain Growth in Metals
 | Surrey University Press, London, 1976
 
 | Humphreys, F. J., Hatherly, M.
 | Recrystallization and Related Annealing Phenomena
 | Pergamon Press, 2003
 | ISBN: 978-0-08-044164-1
 
 | Gottstein, G.
 | Physical Foundations of Materials Science
 | Springer, Berlin, 2010
 | http://dx.doi.org/doi:10.1007/978-3-662-09291-0
 
 | Gottstein, G., Shvindlerman, L. S.
 | Grain Boundary Migration in Metals: Thermodynamics, Kinetics, Applications
 | CRC Press, Boca Raton, 2010
 | ISBN 9781420054354
 
 | Hallberg, H.
 | Approaches to Modeling of Recrystallization
 | Metals, 2011, 1, 16-48
 | http://dx.doi.org/doi:10.3390/met1010016