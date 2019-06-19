# Sweep sim code

This matlab code is for simulating the MRI signal produced in gradient echo MRI using an EPG simulation framework. This code is described in the following papers, if you find this code useful please reference this work. 

> Jackson, LH, Price, AN, Hutter, J et al. Respiration resolved imaging with continuous stable state 2D acquisition using linear frequency SWEEP. Magn Reson Med. 2019; 00: 1– 15. https://doi.org/10.1002/mrm.27834 

> Malik, S. J., Teixeira, R. P. and Hajnal, J. V. Extended phase graph formalism for systems with magnetization transfer and exchange. Magn. Reson. Med. 2018; 80: 767-779. https://doi.org/10.1002/mrm.27040

## Getting Started

Clone this repo onto your local machine and open 'sweep_sim_main'. 

### Prerequisites

This code uses the EPG-X simulation code described by Malik et al. (2018). A barebones version is included wtihin this repo and full code version can be found [HERE](https://github.com/mriphysics)

OPTIONAL: These simulations can be quite intensive, especially when flow is included. The script uses the parfor function to  parallelise much of the simulation. If you are able to connect to a networked machine with more processing power then you can enable the "offload=1" option in sweep_sim_EPG_2.m to greatly accelerate the simulation. This requires the send2remote function [(available here)](https://github.com/laurencejackson/send2remote). 

### Installing

Clone the sweep-sim repo onto your local machine

```
git clone https://github.com/laurencejackson/sweep-sim.git
```

## Running the tests

open sweep_sim_main.m

From here simulations parameters can be tuned and the sims can be run, additional options can be found in sweep_sim_EPG_2.m although these should not need to be changed for most users. 

Simulation parameters are stored in 3 matlab structures

* Tissue - imaging object params
* RF - sequence paramss
* motion - motion params

Code pipleine overview
1. Define simulation paramters
2. Calculate pulse profile
3. Create flipmat - matrix represting applied flip angle with time and space position
4. calculate extended phase graph for every spatial position in flipmat
5. return matrix of spatial signal vs time - s0

## Authors

* **Laurence Jackson** 
* **Shaihan Malik** 
* **Rui Pedro AG Teixeira**
* **Joseph Hajnal**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

