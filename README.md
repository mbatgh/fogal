
# fogal

The code provided here can be used to optimize Force field (FF) parameters for small organic molecules
by fitting them to reproduce energies and forces from ab-initio or DFT calculations.
The resulting FF can be used in Molecular Dynamics simulations. If the input (initial parameters)
are taken from GAFF, and the partial charges are generated accordingly (AM1BCC or RESP), we expect
the resulting FF to be compatible with the Amber family of protein/peptide FFs.
The input  includes, for a given molecule: a topology file with initial values for all parameters, 
a file with coordinates of a series of conformers, and a file that lists energies and/or forces 
corresponding to each conformer. The output is another topology file with optimized force field 
parameters, and a table comparing energies and forces obtained with the initial vs the optimized
FF parameters.

## Features

- Optimization is done using a genetic algorithm, and all the relevant hyper parameters can be adjusted by the user.
- The optimized energy function, and corresponding parameters, are those used in the popular GAFF FF.
- Optimized parameters include: force constants bond-angle, and dihedral force constants, equilibrium
  bond angles, and/or Lennard-Jones paramters epsilon/sigma.
- Input conformers can be a set of single molecules, dimers, trimes, etc.
- Both energies and/or forces can be used for parameter fitting.
- All energy/force values can be assigned individual weights.
- In/output in Gromacs topology file format.

## Prerequisites

The program fogal is a plain C code that requires only standard libraries, and should compile
with gcc on most Linux distributions. Software for the generation of the training data (accurate
energies and forces for the conformers) is not included here, but there are a good number of programs
out there that can be used for this purpose, and some examples are listed below.
Starting from such data the generation and pre-processing of input data is straight forward, but 
execution of the entire pipeline can be facilitated/accelerated with a number of scripts provided
here in the scripts folder, and a number of external tools, including:

- [openbabel](http://openbabel.org/)
- [Gromacs](http://www.gromacs.org/)
- [acpype](https://github.com/alanwilter/acpype)
- [LOOS](https://github.com/GrossfieldLab/loos)
- [Ambertools](http://ambermd.org/AmberTools.php)

For more details see the documented [workflow](examples/dps.md) in the examples folder.

A small, biased, selection of software that can be used for the calculations of training data
(energies and forces for different conformers):

- [Orca](https://orcaforum.kofo.mpg.de/app.php/portal)
- [Psi](https://psicode.org/)
- [Quantum-Espresso](https://www.quantum-espresso.org/)
- [NWChem](https://nwchemgit.github.io/)
- [Gamess-US](https://www.msg.chem.iastate.edu/gamess/)
- [CP2K](https://www.cp2k.org/)

## Installation

```
git clone https://github.com/mbatgh/fogal.git
cd fogal/src
make 
make install
```
should generate a binary called fogal in the bin folder in your home directory.

## Usage

fogal is executed on the command line. Various settings can be controlled in a parameter input file,
or with command line parameters. For a concise overview consider the output of the command fogal -h,
as shown below.

fogal -h

```
# FOGAL version 21



usage: fogal [command-line parameters]

  -h          show this msg
  -z  string  parameter filename
  -o  string  type of parameters optimized [a]ngles,[d]ihedrals,[i]mpropers)
  -g  int     total number of generations
  -i  int     total number of individuals per generation
  -a  float   mutation rate applied to first generation population
  -m  float   mutation rate at the level of a single genome
  -u  float   mutation rate at the level of population
  -c  float   crossover probability (0-1)
  -f  float   fraction of individuals that always get a chance to mate
  -l  float   fraction of poor individuals that get a chance to mate
  -s  int     seed for random generator
  -w  float   cut-off radius for VdW interactions
  -p  string  name of input topology file
  -e  string  name of text file with QM energies
  -t  string  name of trajectory file with conformations
  -k  int     number of conformers
  -q  int     number of QM energies in efw file
  -n  int     number of molecules per conformer
  -v  int     number of forces in efw file
  -r  string  run id, for unique filenames
  -y          test: parameters are read and processed but optimization is not done

```

For more details on the generation of input files, their syntax, and the overall workflow, see the 
example provided [here](examples/dps.md)

## Reference

If you use this code in published work, please cite the reference given below.

Michael Brunsteiner, Sten Nilsoon-Lill, Lucy M. Morgan, Adrian Davis, and Amrit Paudel
Finite temperature mechanical properties of organic molecular solids from classical molecular simulation.
To be submitted (2022).

## Acknowledgements

The work involved in the design of this tool was done at the [Research Center Pharamaceutical Engineering](http://www.rcpe.at),
in a collaborative K project, towards understanding chemical stability of small molecule drugs in the solid state, with funding
from the [Austrian COMET program](https://www.ffg.at/en/comet/programme), and a number of companies (Astra-Zeneca, Janssen,
Pfizer, and UCB).

