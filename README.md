# HVDC

Solves the charge distribution problem on Octree mesh.

The system to be solved is the following:
```math
\begin{cases}
    \frac{\partial \rho}{\partial t}&-\nabla\cdot(\sigma\nabla\varphi)&=0\\
    &-\nabla\cdot(\epsilon_0\epsilon_\infty\nabla\varphi)-\rho+\sum_k \pi_k&=0\\
    \tau_k\frac{\partial \pi_k}{\partial t}&+\nabla\cdot(\epsilon_0\chi_k\nabla\varphi)+\pi_k&=0
\end{cases} 
```
on a cubic domain, for k = 1,2,3.

## Basic usage instructions
### Brief overview
The application is basically composed of a main file `HVDC_main.cpp` in `HVDC/src`, some header files in `HVDC/include`, and some `cpp` files implementing plugins and related necessary structures in `HVDC/src/plugins`. Finally, to properly run the application a `json` parameter file is needed.

The file `tests.h` defines a set of tests describing the problem to solve. The user can add his own test to the file or write it in a new plugin.

Similarly the user can modify the function in the `voltages.h` file or add a new one.

The main program will read from the parameter file
- an output folder where to save results.
- a vector of tests to run.
- the parameters needed by:
  - the main program.
  - the selected tests.

### Installation
The direct dependency of the code is `bim++`, while `Paraview` and `GNU Octave` are necessary
for the post-processing operations. 

For the compilation, a `local_settings.mk` file must be created in order to add or modify, 
accordingly to the user necessities, the compilation options present inside `Makefile`, 
and then the following command must be executed 
```
make
```

### Parameters setup
The application needs a `.json` object file in order to run. Changing variable names inside will cause an exception.
To generate a default parameter file run the command
```
mpirun HVDC_main --generate-params <file-name>
```
If no \<file-name\> is provided the application will generate a file 'data.json'.
The main variables in the file are
```
output_location -> directory where to save results
test_to_run -> vector with names of tests to run in sequence
test1 {} -> object with parameters of test1
test2 {} -> object with parameters of test2
...
```
Each test object has three sections inside:
```
physics
    epsilon_0 -> dielectric constant in vacuum
    physics_plugin -> name of the plugin containing the test
    plugin_params {} -> object containing physical constants used by the test functions
algorithm
    possibly some parameters for the mesh generation like (just as an example)
        NUM_REFINEMENTS -> used to generate a uniform mesh
        maxlevel -> max local level of refinement for a mesh
    T -> Full time: the application will simulate the evolution of the problem from time = 0 to time = T
    initial_dt_for_adaptive_time_step -> starting dt with wich to start the adaptive time step (small value suggested)
    tol_of_adaptive_time_step -> tol of relative error of displacement current
    voltage_plugin -> name of the plugin containing the voltage function
    voltage_name -> name of the function V(t)
    voltage_plugins_params -> contains parameters of the voltage function
    start_from_solution -> if true the application will resume a previous simulation from a temporary solution
    save_temp_solution -> if true the application will save the current solution and allow the user to resume the simulation later
    temp_sol {} -> contains:
        file_of_starting_solution -> (if start_from_solution = true) file name of the solution to resume the simulations
        save_every_n_steps -> to tell how frequently should the application save the temporary solution
options
    save_sol -> if true, the application saves the solution to be exported with octave
    print_solution_every_n_seconds -> time step magnitude; how often the application outputs the solution
    compute_charges_on_border -> if true, the application outputs in a file the different charges on the contact (1)
    save_displacement_currents -> if true, enables the computation of displacement current
    compute_2_contacts -> if true compute displacement current on both the contacts
    save_error_and_comp_time -> if true, tracks error of the adaptive time step and computational times
```


### Command line options

To run the application in parallel with MPI, use the command

```
mpirun -np <number of processes> HVDC_main -f <parameters_file_name.json>
```
If no parameter file is given the application will look for a `data.json` if this exists, otherwise it terminates.

If the application is starting a new simulation (start_from_solution = true) and finds that the directory `<output-folder>/<test-name>` already exists will give a warning before overwriting existing files. To disable this warning use the option `--overwrite` in the command line.

## Post-processing of Output Files
### Visualize the solution

The application output is composed of a series of frames to export and visualize in `Paraview`, and three files containing additional information.

### Visualize the solution

The problem solution at each time step is located in 
```
<output_folder_name>/<test_name>/sol
```
where \<output_folder_name\> is set in the data file.

The user can use the `GNU Octave` function included in `export_phi_rho_p1_3.m` 
(located in the  folder `script/m`) to generate the .vtu files,

which can by opened using Paraview. This function is a wrapper over the function `export_tmesh_data.m`
provided with bim++, so it is necessary to add the path `script/m` of bim++ to the Octave path. 
The usage is the following:
From the terminal:
```
cd <output_folder>/<test_name>/sol
octave
```
In octave run
```
export_phi_rho(<initial time step>, <final time step>, <number of processes>)
```
with `<initial time step>` and `<final time step>` being the indexes of the first and last time step and `<number of processes>` 
equal to the number of processors. 

### Other files

If the correspondent options are enabled the application ouptuts the files
```
charges_file.txt
I_displ_file.txt
error_and_comp_time.txt
```
in the folder `<output_folder>/<test_name>`.

By default the application will export these files also in `.json` format.

It is possible anyway to explicitly export a (txt) file in `json` format by the command
```
mpirun HVDC_main --export-json <file/name.txt>
```
To import the data in `Octave` we use the script in `script/m/import_json_file.m`.
