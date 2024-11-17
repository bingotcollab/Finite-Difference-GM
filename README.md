The Finite-Elements-GM is a code dedicated to calculate the Quasi-Normal Modes (QNM) of a scalar field over a Gutsunaev-Manko (GM) background.
The Code consists of two main classes:

* **Metric**: Defines the GM metric, i.e., its parameters, as well as the functions g and f, defining the coefficients in the Weyl gauge.
* **Solver**: Defines the functions required for solving the field equations. It receives the metric and the simmulational parameters. It cotains the method *solve*,
which iterates the finite-elements method to evolve the field over time, and periodically saves the results. If the load parameter is passed, it will load the last result
from the previous execution, and start the iteration from there.

For further details, see E.C. Ribeiro, L. Formigari, M.R. Ribeiro Jr., E. Abdalla, et all, Stability of the spacetime of a magnetized compact object.  [[arXiv:2411.nnnnnn]](https://arxiv.org/abs/2411.nnnnn).


**Execution**

To run the code, the argument --params should be passed followed by the path to the parameters json file. If loading from a previous run, pass the argument --load.

Example:

```
  python finite-differece-gm.py --params ./example_parameters.json
```

**Parameters**

   The parameters file should pass all values needed to execute the code. The example_parameters.json contain all the ones needed to run the code as is in this repository. The meaning of the parameters are as follows.

  * metric: A dictionary containing all parameters that deffine the used metric. For GM, there's a single one "alpha", the magnetization constant.
  * lmax: Number of degrees of the spherical harmonics to be used in the simulation
  * m: The order m of the spherical hamonics to be used.
  * Npoints_x: Number of points in the grid for the radial coordinate 
  * Npoints_t: Number of points in the grid for the time coordinate
  * xmin: Position of the lower bound of the grid in the radial coordinate. For the code as is, it corresponds to a mirror at this location.
  * xmax: Position of the upper bound of the grid in the radial coordinate. For the code as is, it corresponds to a mirror at this location.
  * t_step: Time interval between consecutive steps of the simulation.
  * initial_distance: Position of the peak of the initial gaussian excitation.
  * initial_spread: Twice the variance of the initial gaussian excitation.
  * li: Degree of the spherical harmonic in which the initial excitation occours.
  * observing_points: Number of points to be sampled and saved at the outputs.
  * block_size: Number of iterations done before an intermediate saving step. The solver will allocate a matrix of size lmax\*Npoints_x\*block_size, so reducing this parameter reduces the RAM used, at the cost of an slower execution.
  * save_full: If selected, all simulated points will be saved to output files.
  * charge: charge of the scalar field.

**Outputs**

  The code generates three output files:

  * Psi_from_t_SSSS_to_EEEE.npz: File containing the values of the field at the observation points from the time SSSS until the time EEEE, compressed to an .npz file.
  * last_two_times_TTTT.npz: Contains the value of the field in all positions for two consecutive times, also compressed to an .npz file. It is used by the solver to reload a previous simulation if the parameter --load is passed.
  * observation_points.npz: Contains the radial positions of the observation points, in turtle coordinates, also compressed to an .npz file.


