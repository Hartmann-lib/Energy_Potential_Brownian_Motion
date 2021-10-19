# Energy Potential Diffusion Simulation

In literature, the heterogeneous conformations of a biomolecule (e.g. protein or DNA) are often represented with an one-dimensional energy landscape. In this landscape conformational states are energy wells separated by barriers. Interconversion dynamics between these states are then described by conformational diffusion within the energy landscape. 

In this simulation package, several types of one-dimensional energy landscapes are provided such as an harmonic, quartic, sinusoidal and custom built potential. State trajectories are simulated by Langevin dynamics using the gradient of the selected energy landscape, grad(*G*(*x*)), and the diffusion constant, *D*. The script also evaluates the passage time between two reaction coordinates, enabling for example, the detection of the transition path time, which denotes the time needed to overcome the energy barrier.

<p align="center">
   <img src="https://user-images.githubusercontent.com/58071484/137738308-d6881a79-d1a1-4096-99b3-b96e57e12e9e.JPG" height=50>
</p>

### Selection of the energy landscape

Choose between an harmonic, quartic, sinusoidal or custom built energy landscape by selecting the corrsponding model function in the main script. Feed the model function with the necessary features of the energy landscape.

![simNRJDiff_Figure1](https://user-images.githubusercontent.com/58071484/137720336-c499caca-533c-4e00-8c06-2379855c89da.png)

### Simulation output
After performing the simulation, the reaction coordinate histogram (grey bars, left panel) collected from all trajectories is plotted. The similarity of the histogram with the probability density function, *p*(*x*), derived from the energy landscape (black line), shows if the landscape is sufficently sampled, hence &#916;*t* is small enough. In the right panel the passage time histogram between *q*<sub>L</sub> and *q*<sub>R</sub> is plotted. The bottom panel depicts the last state trajectory. All generated trajectories can be saved in a directory by setting *boolSave*=1 (see example files in folder "data\\").

<p align="center">
  <img src="https://user-images.githubusercontent.com/58071484/137707404-58e4e83a-afaf-4015-bbdd-ccb9bb040450.png">
</p> 

**Note:** In the script above I make use of the optimization function *fminsearchbnd* to generate a custom energy landscape. The author of *fminsearchbnd* is John D'Errico, who released the function 2006 on FileExchange - Mathworks (https://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon?s_tid=srchtitle).
