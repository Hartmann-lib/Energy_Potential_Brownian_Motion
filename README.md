# Energy Potential Diffusion Simulation

In literature, the heterogeneous conformations of a biomolecule (e.g. protein or DNA) are often combined in an one-dimensional energy landscape. In this landscape conformational states are wells separated by energy barriers. Interconversion dynamics between these conformational states are then described by diffusion within the energy landscape. 

In this simulation package several types of one-dimensional energy landscapes are provided such as a harmonic, quartic, sinusoidal and custom built potential. State trajectories are simulated by Langevin dynamics using the gradient of the selected energy landscape. The script also evaluates the passage time between two reaction coordinates, enabling e.g. the detection of the transition path time, the time needed to overcome the energy barrier.

### Selection of the energy landscape

Choose between an harmonic, quartic, sinusoidal or custom built energy landscape by selecting the corrsponding model function in the main script. Feed the model function with the necessary features of the energy landscape.

![simNRJDiff_Figure1](https://user-images.githubusercontent.com/58071484/137707374-0c9cbd84-a050-4cb5-a97e-8b089bccb3c6.png)


![simNRJDiff_Figure2](https://user-images.githubusercontent.com/58071484/137707404-58e4e83a-afaf-4015-bbdd-ccb9bb040450.png)
