The .cpp and .h files can be complied with the LAMMPS package to implement a classical force field for vibrationally excited water molecules, based on the spectroscopic maps of Skinner and co-workers.

The derivation of the force field is documented in mapForce.pdf.

To use the force field, copy the .cpp and .h files to the <lammps>/src directory and compile LAMMPS.
The force field defines a new pair_style in LAMMPS: ```pair_style excited/map```.
When the pair_style is specified in the LAMMPS input script, it takes five arguments

```pair_style 	excited/map <cutoff> <oxygen type> <special oxygen id> <mapA> <mapB>```

- ```cutoff``` is the cutoff for computing the electric field for the map (distance units)
- ```oxygen type``` is the LAMMPS type that corresponds to oxygen
- ```special oxygen id``` The LAMMPS id of the oxygen atom on the molecule that will have the excited stretch. The next hydrogen atom in order (id+1) is chosen as the excited hydrogen atom.
- ```mapA``` and ```mapB``` are the parameters of the map. The map has the form omega=mapA\*E + mapB\*E^2 + constant (see mapForce.pdf for details).
     - ```mapA``` has LAMMPS units of charge\*length.
     - ```mapB``` has lammps units of (charge\*length)^2/energy.

The pair style is meant to be used in conjunction with pair_style hybrid/overlay, but this is not strictly necessary.
Several example input scripts illustrating this are given in the inputScripts/ directory. These scripts use the rigid SPC/E water model. The three relevant scripts are
- in.gs: water in the vibrational ground state (i.e. normal SPC/E water)
- in.es: water with one decoupled OH stretch in the excited state. This implicitly treats all the water molecules as D2O molecules, with the one special one treated as an HOD molecule.
- in.es_neq: the system is equilibrated in the vibrational ground state, then the OH stretch is excited and the system relaxes to this perturbation.

The other input scripts are used within those three and do not stand alone.
