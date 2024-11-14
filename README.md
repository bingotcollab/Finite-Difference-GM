The Finite-Difference-GM is a code dedicated to calculate the Quasi-Normal Modes (QNM) of a scalar field over a Gutsunaev-Manko (GM) background.
The Code consists of two main classes:

* **Metric**: Defines the GM metric, i.e., its parameters, as well as the functions g and f, defining the coefficients in the Weyl gauge.
* **Solver**: Defines the functions required for solving the field equations. It receives the metric and the simmulational parameters. It cotains the method *solve*,
which iterates the finite-elements method to evolve the field over time, and periodically saves the results. If the load parameter is passed, it will load the last result
from the previous execution, and start the iteration from there.

The Gutsunaev-Manko Metric is stated in the Weyl Gauge, with the following line element:


**Execution**
To 


**Outputs**



**Further References**
