How to get? How to compile?
===========================

1. Check your Linux installation for a working installation of **cmake** and **make**
2. Pull the project from the git repository https://github.com/mkuehbach/SCORE
3. Check your compiler. Both the Intel compiler (v14.0 tested) and the GNU compiler (v4.8.3 tested) can be utilized.
4. Make sure there is a folder with a **src**, a **build** folder, and the CMakeLists.txt file.
5. Make sure that the source code files are in **src** while a **SCORE.*.uds** parameter file is in **build**.
6. Open a console and type *cd build*
7. Only once when setting up a new computer *cmake -DCMAKE_BUILD_TYPE=Release ..* (up to O2 optimization tested)
8. Compile the program with *make*
9. Find happily the binary in the **build** folder