# This submit a job in bigred2

JOBSCR = macroassembly.pbs
TESTSCR = test.pbs
E=-1500
K=-600
r=1
n=100 
c=0.150 
q=45 
D=56 
d=6.7

e=50
t=1
x=1
i=0 
N=2001

NODESIZE=4


run-preprocessor:
	./preprocessor -E $(E) -K $(K) -n $(n) -c $(c) -r $(r) -q $(q) -D $(D) -d $(d)

run-postprocessor:
	./postprocessor -e $(e) -t $(t) -x $(x) -i $(i) -N $(N)

submit: 
	@echo "Launching the job";
	qsub $(JOBSCR)

run-local-parallel:
	time mpirun -n $(NODESIZE) lmp_g++ < in.lammps

run-local-serial:
	time lmp_g++ < in.lammps

