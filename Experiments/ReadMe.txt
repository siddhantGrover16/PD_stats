The Experimental evaluation of our project is presented in this directory

1.  GenerateTrees.py uses dendrpy to create birthdeath trees
2*. randvsavg_0.1.py compares the runtime of the avg algorithm to randomly sampling taxa till a 0.1% convergence is observed
3*. randvsavg_0.5.py compares the runtime of the avg algorithm to randomly sampling taxa till a 0.5% convergence is observed
4*. randvsvar_0.5.py compares the runtime of the var algorithm to randomly sampling taxa till a 0.5% convergence is observed
5.  The Scalability_study directory presents a basic scalability study for trees ranging from sixe 1000 to 10000 where k=n. The
   trees in the section are created via the dendropy method mentioned in GenerateTrees.py

* this program uses scikit-bioto compute the pd of sampled sets. The experiments are run on t500.tre amd t500_.tre(also present in this directory for testing)

