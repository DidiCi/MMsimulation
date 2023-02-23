# MMsimulation
### Description
A Monte Carlo method coupled to a genetic optimization algorithm to simulate DNA replication and compare the results to experimental data.
The method was first described and used in the following study: [Ciardo et al. 2021](https://www.mdpi.com/2073-4425/12/8/1224).
The model MM5, described in the article, is used to simulate DNA replication.
### Installation
To use the code, download the repository on your computer and unpack it. Open the repository in Matlab as "Current Folder".
### Hardware and OS Requirements
The simulation and fit require a standard computer with enough RAM to support the in-memory operations. 
The scripts have been tested on a computer with the following specifications: 6 cores, 2.60 GHz, 32GB RAM.
The scripts have been tested on the following systems: Windows 10 and Ubuntu 20.04.
### Software dependencies
The scripts have been tested with the following Matlab versions: R2020b.

The following Matlab toolboxes are required:
```
Global Optimization Toolbox
System Identification Toolbox
Statistics and Machine Learning Toolbox
```

 
# Demo
In order to test the scripts, a small dataset is provided in the folder "Data_demo/data_example". In this folder, each file contains the position and raw intensity measurements for a single combed DNA molecule. The file "Log.txt" contains the list of file names. 
The scripts must be executed in the following order:
1. `Data_extraction/fittotot.m`: this script is used to extract, process and save the experimental data in a structure array to be used later;
2. `Simulation/geneticalgorithm_main.m`: this script simulates DNA replication and compares results with experimental data to optimize simulation parameters; the script can perform multiple rounds of optimization;
3. `Result_analysis/finalanalysis.m`: this script selects for each optimization round the best individual and calculates the parameters to be used for images and statistics;
4. `Images/analysis1`: this script plots the following parameters for experimental and simulated data: replicated fraction f(t), rate of origin firing and fork density were calculated for each molecule as a function of f(t), eye-to-eye distances, eye and gap length distributions; the statistical tests are printed in an excel file;
5. `Images/analysis2`: when multiple conditions must be compared, this script plots the mean and standard deviations of the simulation parameters for the two conditions; the statistical tests are printed in an excel file.

The second script is the most time-consuming. The expected execution time for an optimization round is about 45min by using 4 workers in parallel and a computer with the specifications listed above.
