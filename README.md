# MATLAB Codes for Variable-Order Nonlocal Equations

This repository contains the MATLAB codes used to reproduce the numerical experiments in the manuscript:

**“Optimal error estimates of a collocation method for variable-order nonlocal equations with graded meshes.”**

## Repository Structure

The four folders correspond to different numerical experiments:

| Folder                                              | Description                                                                                                                                                                                             |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`Constant_Frac_Lap_vpa`](./Constant_Frac_Lap_vpa/) | Codes for the constant-order fractional Laplacian benchmark problem. Variable-precision arithmetic (VPA) is used to construct strongly graded meshes and accurately evaluate the discrete coefficients. |
| [`variable_Lap_ref`](./variable_Lap_ref/)           | Codes for variable-order nonlocal problems in which the exact solution is unavailable. A numerical solution computed on a sufficiently fine mesh is used as the reference solution.                     |
| [`variable_Lpa_MMS`](./variable_Lap_MMS/)           | Codes for the variable-order examples based on the method of manufactured solutions (MMS). The source term is evaluated numerically from the prescribed manufactured solution.                          |
| [`variable_mixed`](./variable_mixed/)               | Codes for the numerical experiments involving mixed local and nonlocal operators. These experiments illustrate the applicability of the proposed collocation method to a broader class of models.       |

Each folder contains the MATLAB functions and demonstration scripts required for the corresponding experiment.

## Requirements

* MATLAB
* Symbolic Math Toolbox, required for the VPA computations in `Constant_Frac_Lap_vpa`

## Running the Codes

1. Download or clone this repository.
2. Open MATLAB.
3. Set the selected experiment folder as the current MATLAB working directory.
4. Run the corresponding demonstration script in that folder.
5. The maximum errors and convergence rates will be displayed in the MATLAB Command Window.

For computations on strongly graded meshes, especially when the fractional order is small, the VPA-based implementation should be used to prevent the coalescence of mesh points caused by double-precision rounding.

## Notes

* Both uniform and graded meshes are considered in the numerical experiments.
* The grading parameter is selected according to the setting described in the manuscript.
* Fine-mesh reference solutions are used when closed-form exact solutions are unavailable.
* Computation time may increase substantially for fine meshes and VPA-based calculations.
