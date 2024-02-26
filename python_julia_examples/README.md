In this folder, we show two examples on how to use Julia TN-AHS simulator in combination with Braket Python SDK:

* `introduction.ipynb`: Introduction notebook [this Braket example](https://github.com/amazon-braket/amazon-braket-examples/blob/main/examples/analog_hamiltonian_simulation/00_Introduction_of_Analog_Hamiltonian_Simulation_with_Rydberg_Atoms.ipynb) using the tensor network simulator. We have added a section at the end of the notebook to show how to save a Braket analog Hamiltonian simulation (AHS) program into a JSON file, followed by using python subprocess to run the program with Julia-based tensor network simulator. 

* `mis_problem.ipynb`: Solving Maximum Independent Set (MIS) problem on a square nearest neighour lattice using piecewise linear annealing protocol.


The final results are saved in the folder `experiments/`. 

A useful guide for how to install Julia may be found at [here](https://ferrolho.github.io/julia/linux/ubuntu/how-to-install-julia-on-ubuntu/).