# MMsimulation
### Description
A Monte Carlo method coupled to a genetic optimization algorithm to simulate DNA replication and compare the results to experimental data.

### Installation
To use the code, download the repository on your computer and unpack it. Open the repository in Matlab as "Current Folder".
### Hardware and OS Requirements
The simulation and fit require a standard computer with enough RAM to support the in-memory operations. 
The code has been tested on a computer with the following specifications: 6 cores, 2.60 GHz, 32GB RAM.
The code has been tested on the following systems: Windows 10 and Ubuntu 20.04.
### Software dependencies
The following Matlab functions/packages are required:
```
ciao
```

 
# Demo
In order to test the scripts, a small dataset is provided in the folder "Data_demo/data_example". In this folder, each file contains the position and raw intensity measurements for a single combed DNA molecule. The file "Log.txt" contains the list of file names. 
The scripts must be executed in the following order:
1. `Data_extraction/fittotot.m`:
2. `Simulation/geneticalgorithm_main.m`:
3. `Result_analysis/finalanalysis.m`:
4. `Images/analysis1`:
5. `Images/analysis2`:

Expected about 45min for one optimization of the simulation with 4 cores in parallel.
