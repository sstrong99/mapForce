#LAMMPS input file
#to simulate bulk SPC/E water

#####################################################################
clear

variable	samp_rate 	equal 1
variable	thermo_rate 	equal 1
variable	equil		equal 1000
variable	run		equal 1000

variable	ts		equal 2.0
variable	Tdamp		equal 100*${ts}
variable	myT		equal 298.0
variable	Wlat		equal 3.1074  #for water density 0.997

units 		real
atom_style	full

dimension 	3

boundary 	p p p

#############################################################################
#setup box

lattice	  	sc ${Wlat}
region		simbox block -3 3 -3 3 -3 3 units lattice

#############################################################################
#set up potential

create_box	2 simbox bond/types 1 angle/types 1 extra/bond/per/atom 2 &
		extra/special/per/atom 2 extra/angle/per/atom 1

molecule 	spce spce.mol
mass            1 15.9994 #oxygen
mass            2 1.008   #hydrogen

pair_style 	lj/cut/coul/long 9.0
pair_modify 	table 0 table/disp 0 shift yes

bond_style  	harmonic
angle_style 	harmonic

kspace_style 	pppm 1.0e-6

pair_coeff 	* * 0.0 0.0
pair_coeff	1 1 0.1553 3.166

#potential coeffs very stiff for minimization
bond_coeff	1 100000 1.0
angle_coeff	1 100000 109.47

#############################################################################
#setup for run
thermo 		${thermo_rate}
thermo_style 	custom step vol temp epair pe etotal press density

timestep 	${ts}
run_style 	verlet

neighbor	2.0 bin
neigh_modify	every 1 delay 3 check yes

#############################################################################
#make atoms and rough minimze
create_atoms	0 box mol spce RAND

min_style 	cg
minimize	1.0e-4 1.0e-4 10000 100000

#potential coeffs aren't too important since will be rigid anyways
bond_coeff	1 554.13 1.0
angle_coeff	1 45.769 109.47

#############################################################################
#initialize velocity and rigid constraint

fix 		rigid all shake 1.0e-8 100 0 b 1 a 1 t 1 2 mol spce
velocity 	all create ${myT} RAND dist gaussian rot yes mom yes

#scale velocity
run 0
velocity all scale ${myT}

#############################################################################
#equilibrate bulk water with NVT

fix 		1 all nvt temp ${myT} ${myT} ${Tdamp}
fix_modify	1 energy yes

run 		${equil}
velocity all scale ${myT}

#############################################################################
#run at NVE

unfix           1
fix 		2 all nve

dump 		1 all custom ${samp_rate} dump.lammpstrj id x y z type
dump_modify	1 sort id
dump 		2 all xtc ${samp_rate} dump.xtc  #sorts by default

run 		${run}