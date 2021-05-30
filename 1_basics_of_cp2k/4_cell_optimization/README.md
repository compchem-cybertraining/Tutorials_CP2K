# Cell optimization

More explanation has been given in the geometry optimization section. Here, we just mention some notes on cell optimization. It is better to perform both cell and geometry 
optimization simultaneously. This is done by setting the `TYPE` in `&CELL_OPT` to `DIRECT_CELL_OPT`. There are different types of cell optimization such as `MD` which is not considered here. 

As was mentioned in the geometry optimization section, we need accurate forces
so that the optimizer can better optimize the structure. Therefore, more accurate electronic structure calculations are required which is done by lowering the value of `EPS_SCF`. The `STRESS_TENSOR` is also needed to compute the pressure for cell optimization. You can also print out the stress as well
using `&STRESS ON` in `&PRINT` section.

Here, we did not perform cell optimization because the changes in the geometry optimizations were small. But you can do it if the structure in the geometry optimization
shows massive changes.
