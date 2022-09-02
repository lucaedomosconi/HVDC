# HVDC

Solves the charge distribution problem on Octree mesh.

Solves the following system of PDEs 

\f{equation}{ 
\begin{cases}
    -\nabla\cdot(\epsilon_0\epsilon_\infty\nabla\phi) =\rho\\
    \frac{\partial \rho}{\partial t}+\nabla\cdot(\sigma\nabla\phi)=0
\end{cases} 
\f}

on a cubic domain.

## Basic usage instructions

### Installation
The direct dependency of the code is `bim++`, while `Paraview` and `GNU Octave` are necessary
for the post-processing operations. 

For the compilation, a `local_settings.mk` file must be created in order to add or modify, 
accordingly to the user necessities, the compilation options present inside `Makefile`, 
and then the following command must be executed 

```
make all
```

### command line options

To run the code in serial, the necessary commands are

```
mkdir results
cd results
../HVDC_main 
```
 
To run in parallel with MPI, change the last line with

```
mpirun -np <number of processes> ../HVDC_main
```

with `<number of processes>` equal to the desired number of processors. 

To select a different test case, select the appropriate header file to include in HVDC_main.cpp.
The available files are:

```
test_1.h -> uniform epsilon in the domain
test_2.h -> two values of epsilon for each half of the domain
test_3.h -> three paper layers and two oil layers with oil filled holes
test_4.h -> three paper layers and two oil layers with air filled holes
test_5.h -> spherical air filled hole
```

### Post-processing of Output Files

The user can use the `GNU Octave` function included in `export.m` 
(it is located in the  folder `script/m`) to generate the .vtu files,
which can by opened using Paraview. The usage is the following:
```
export(<number of time steps>,<number of processors>)
```
Moreover,  the file `rho.m` contains a function that generates a plot
of the value of the charge density rho over time. The usage is the following:
```
rho(<refinement level>,<cube side lenght>)
```
where the refinement level refers to the elements selected by the `rho_idx` variable,
and the cube side length is usually set to 0.001.
