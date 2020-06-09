# This submit a job in bigred2

JOBSCR = macroassembly.pbs
TESTSCR = test.pbs
Q=1500
n=25 
c=0.150 
q=35 
D=56 
d=6.7
i=0 
N=150
NODESIZE=4

run-preprocessor:
	./preprocessor -Q -$(Q) -n $(n) -c $(c) -q $(q) -D $(D) -d $(d) -i $(i) -N $(N)

submit: 
	@echo "Launching the job";
	qsub $(JOBSCR)

run-local-parallel:
	time mpirun -n $(NODESIZE) lmp_g++ < in.lammps

run-local-serial:
	time lmp_g++ < in.lammps

