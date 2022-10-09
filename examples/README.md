
Below an example workflow for an optimization of FF parameters for diphenyl-sulfone is provided

1. Starting with a structure file, e.g. an sdf file or a smiles string, generate a conformer of the molecule in xyz format with canonical atom ordering. Using openbabel, we'd say:

```
obabel dps.sdf -Oc-dps.xyz --canonical
```

2. Generate multiple conformers, for example by using some semi-empirical code like MOPAC or XTB (which we recommend) to perform an MD simulation in vacuum, saving structures at regular intervals. After optimizing the structure (ideally at the same level of theory that is used to calculate enegies and forces) generate a single input structure (here called b973copt-dps.xyz), and use xtb as in:

```
xtb b973copt-dps.xyz --input mdxtb.inp --md --namespace dps > dps.out
```

The resulting trajectory is then converted into individual xyz files, using openbabel, or a simple awk script, to generate individual structure files in xyz format, as in:

```
awk '{if(NF==1) {na=$1; nc++; fn=sprintf("dps-%05d.xyz",nc);} print >> fn; }' dps.xtb.trj
```

providing:

```
dps-00001.xyz, dps-00002.xyz, etc...
```

3. Optionally cluster the generated conformers to reduce the total number of structures. The LOOS library is one of many tools that can be used for this purpose. For a typical small organic molecule with a limited number of rotatable bonds a few thousand conformers should do for a start, but validating the results by comparing parameter values as obtained when increasing the numbers of conformers is certainly a good idea.

4. Use the generated conformers as input for some program to calculate reasonably accurate energies and forces. Some DFT appoach with a decent basis set and functional usually provides a good prize/performance ratio. Anything beyond double hybrid DFT is probably excessive. We like to use ORCA for this step.

```
orca b973c-dps-0001.inp > b973c-dps-0001.out
```

Which (if you use ORCA) gives you energies and forces on each atom as a file in engrad format.

5. Convert the output of the DFT calculations to a format that fogal can read. Here we use simple awk scripts that can process engrad files as output by ORCA. For other DFT or ab-initio program output you need to adjust the script accordingly (or write you own).

```
./engrad2fogal.awk orca-00000.engrad
mv structures structures-dps.xyz
cat energies forces > efw-dps
```

Now we have a file containing the coordinates of 1024 conformers in xyz format,
and a file contiaining the corresponding 1024 energies, followed by the 1024*25*3 forces
on each atom.


6. To generate a Gromacs topology input file for the molecule with GAFF2 parameters, you'd generate a mol2 type file, and use acpype, as in:

```
obabel dps.xyz -O dps.mol2 
acpype -i dps.mol2 -c user -a gaff2
```

this should create a number of files including:

```
dps_GMX.itp
```

which contains all force field parameters for the molecule in Gromacs itp format.

7. Generate partial charges (these are taken as input an not optimized here). Unless there is a very good reason for not doing so, we recommend using either AM1BCC, RESP, or ABCG2 charges. These can be generated with a number of tools, e.g., the ambertools package. Save the charges in a simple text file, one charge per line, in the same order as in the structure files.

```
dps.chg
```

8. Edit the topology to include the previously calculated charges and account for symmetry constraints. Here we label each force field parameter by an integer number. Parameters that should be identical due to symmetry are assigned the same number. The numbers are included after a semi-colon on each line defining a parameter (or two in the case of bond angles)

```
./print-sym.py dps.pdb > dps.sym
./update_gmx_top_sym.awk -v resname=DPS dps.chg dps.sym dps_GMX.itp b973copt-dps.xyz > sym-abcg2-dps.top
```

9. Finally apply fogal to optimize the parameters of the force field to better reproduce the DFT energies:

```
fogal -z fogal1.inp -o adi -g 50 -i 100 -s 1 -p sym-abcg2-dps.top -e efw-dps -t structures-dps.xyz -n 1 -k 1024 -q 1024 -v 76800 -r fogaldps > dps-stdout 2>&1 &
```

This takes about one minute on a mediocre laptop. For a small highly symmetric molecule like dps this should provide close to converged parameters. For larger/more complex molecules, you obviously need to run tests, varying the total number of conformers, the number of generations (g) and individuals (i) of the genetic algorithm. If you feel adventurous, you can also play with the settings in fogal1.inp which set various other parameters of the genetic algorithm, such as mutation rates, cross-over probability, etc.

10. Visualize results. The output of fogal includes the development of the fittness (normalized difference between DFT and MM energies and forces) as the population evolves. This is written to STDOUT and in the above example redirected to a file dps-stdout which can directly be viewed using, e.g., gnuplot

```
gnuplot> plot 'dps-stdout' u 1:4
```

The energies and forces DFT vs MM are provided in the file fogaldps-ef-ga. For visualization this can be split into energies and forces

```
awk '{if(NR<=1024) print > "e-"FILENAME; else print >> "f-"FILENAME}' fogaldps-ef-ga
```

and, again, with gnuplot:

```
gnuplot> plot 'e-fogaldps-ef-ga' u 5:3
gnuplot> plot 'e-fogaldps-ef-ga' u 5:6
```

This comparsion is shown in the figure below, which demonstrates that with the final parameter set (right) the energy and force values are clearly closer to the DFT based values, than those obtained when using a generic force field like GAFF2 (left).

![optimization](./opt.png)



