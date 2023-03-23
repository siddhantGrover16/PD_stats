from dendropy.simulate import treesim
"""
Dendropy is used to simulate the birth-death trees used for all the experiments done for the project
"""

t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=500)
"A birth death tree with 500 taxa is created and is stored in both newick and nexus format "

t.print_plot()# enables the user to visualize the
t.write(path="t500.tre", schema="nexus") # write tree to file in nexus format
t.write(path = "t500_.tre",schema = "newick") # write tree to file in newick format