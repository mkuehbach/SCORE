Further details for developers.
===============================

In contrast to the pure MPI version in which each automaton domain expands completely in one implicit array (the **mycellgrid**), 
the hybrid version subdivides statically this grid in the y and the z direction. Each such automaton region is now handled by the 
**caregionMemHdl** class which allocates storage and bookkeeps two interface lists. According to this list, cells are either categorized as located completely in the simulation domain or close to its boundary, which is a one cell thick layer.

The cells that are close to this boundary require synchronization with other regions to avoid data races and to maintain a consistent
view of the infection states. For this **halo** regions are defined, into which the threads first infect in parallel, then synchronize **syncHalo**
and again read out in parallel the changes that occurred in the halo regions and require synchronization in the corresponding memory regions.
This minimizes thread synchronization.

At the moment the cell list defragmentation functionality is disabled.