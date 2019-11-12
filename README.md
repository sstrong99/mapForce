The .cpp and .h files can be complied with the LAMMPS package to implement a classical force field for vibrationally excited water molecules, based on the spectroscopic maps of Skinner and co-workers.

The force field is documented in the .tex file.

To use the force field, copy the .cpp and .h files to the <lammps>/src directory and compile LAMMPS.

Several example input scripts are given in the inputScripts/ directory. These scripts use the rigid SPC/E water model. The three relevant scripts are
in.gs: water in the vibrational ground state (i.e. normal SPC/E water)
in.es: water with one decoupled OH stretch in the excited state. This implicitly treats all the water molecules as D2O molecules, with the one special one treated as an HOD molecule.
in.es_neq: the system is equilibrated in the vibrational ground state, then the OH stretch is excited and the system relaxes to this perturbation.

The other input scripts are used within those three and do not stand alone.
