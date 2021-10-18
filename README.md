# Energy Potential Diffusion Simulation

In literature, the heterogeneous conformations of a biomolecule (e.g. protein or DNA) are often represented with an one-dimensional energy landscape. In this landscape conformational states are energy wells separated by barriers. Interconversion dynamics between these states are then described by conformational diffusion within the energy landscape. 

In this simulation package, several types of one-dimensional energy landscapes are provided such as an harmonic, quartic, sinusoidal and custom built potential. State trajectories are simulated by Langevin dynamics using the gradient of the selected energy landscape, grad(*G*(*x*)), and the diffusion constant, *D*. The script also evaluates the passage time between two reaction coordinates, enabling for example, the detection of the transition path time, which denotes the time needed to overcome the energy barrier.

<p align="center">
  <img src="https://user-images.githubusercontent.com/58071484/137738308-d6881a79-d1a1-4096-99b3-b96e57e12e9e.JPG">
</p>

### Selection of the energy landscape

Choose between an harmonic, quartic, sinusoidal or custom built energy landscape by selecting the corrsponding model function in the main script. Feed the model function with the necessary features of the energy landscape.

![simNRJDiff_Figure1](https://user-images.githubusercontent.com/58071484/137720336-c499caca-533c-4e00-8c06-2379855c89da.png)

### Simulation output
In the end of the simulation, the script shows the reaction coordinate histogram (grey bars, left panel) collected from all trajectories. The similarity of the histogram with the probability density function, *p*(*x*), derived from the energy landscape, proves if the landscape is sufficently sampled, hence &Delta*t* is small enough.

the whole energy landscape. In the right panel the passage time histogram between *q*<sub>L</sub> and *q*<sub>R</sub> is shown. The bottom panel depicts the last state trajectory.

<p align="center">
  <img src="https://user-images.githubusercontent.com/58071484/137707404-58e4e83a-afaf-4015-bbdd-ccb9bb040450.png">
</p> 
