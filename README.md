# RAMODO
Reformulation Approaches for Multi-objective Discrete Optimization

This repo is created to include the supporting material for the TUBİTAK project 122M422 - "Reformulation Approaches for Multi-objective Discrete Optimization Problems" and will not be maintained in the future. To solve multiobjective optimization problems, use [MultiObjectiveAlgorithms.jl](https://github.com/jump-dev/MultiObjectiveAlgorithms.jl).

You can load problems and solve them using the algorithms provided.

First, run the main script `MODO.jl`. After that, load the problem using the data structures provided. Code works with the formats given in `txt` folder for each problem set.

```
include("MODO.jl")
```

## Knapsack Problem
$$
\begin{align*}
\max & \sum_{i=1}^p v_{i}^k x_i, \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n w_i x_i \le W \\
& x_i \in \{0, 1\}
\end{align*}
$$

```
kp = MOKnapsackProblem("Kirlik & Sayın/KP/txt/KP_p-3_n-10_ins-1.txt")
```

## Assignment Problem
$$
\begin{align*}
\max & \sum_{i=1}^p c_{ij}^k x_{ij}, & \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n x_{ij} = 1, & \quad \forall j \in 1,...,n \\
& \sum_{j=1}^n x_{ij} = 1, & \quad \forall i \in 1,...,n \\
& x_{ij} \in \{0, 1\}
\end{align*}
$$

```
ap = MOAssignmentProblem("Kirlik & Sayın/AP/txt/AP_p-3_n-5_ins-1.txt")
```

## Integer Linear Problem
$$
\begin{align*}
\max & \sum_{i=1}^p c_{i}^k x_i, & \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n a_{ji} x_i \le b_j, & \quad \forall j \in 1,...,m \\
& x_i \in \{0, 1\}
\end{align*}
$$

```
ilp = MOIntegerLinearProblem("Kirlik & Sayın/ILP/txt/ILP_p-3_m-20_n-10_ins-1.txt)
```

Now, to solve each problem, just pick an algorithm and run it. Note that Kirlik & Sayın also expect an objective funtion index to build two-stage $\varepsilon$-constraint models. Each function also returns statistics about the runs.

Kirlik & Sayın (2014):
```
Y_N, _ = kirlik(kp, 1)
```

Dominguez-Rios (2021):
```
Y_N, _ = dominguez(ap)
```

Tamby & vanderpooten (2021):
```
Y_N, _ = tamby(ilp)
```

`RepresentationMetrics.jl` provides functions for quality measures for representations. They accept the representation and the nondominated set as the arguments. Assuming `Y_R` is the representation somehow obtained, we can compute the following quality metrics:

1. Overall Nondominated Vector Generation Ratio (ONVGR) `ongvr(Y_R, Y_N)`
2. Additive Epsilon Indicator ($\varepsilon_+$) `additive_epsilon_indicator(Y_R, Y_N)`
3. Coverage Error `coverage_error(Y_R, Y_N)`
4. Range Ratio `range_ratio(Y_R, Y_N)`
5. Hypervolume (see the [link](https://lopez-ibanez.eu/hypervolume) for computing $hv(\mathcal{Y}_R)$
`hypervolume(Y_R, reference)`
6. Hypervolume Ratio

In order to use the function to compute hypervolume, you first need to build a shared library (extensions are `.so` for Linux, `.dll` for Windows and `.dylib`for Mac) using the code provided with link.
 
## References
[1] Kirlik, G., & Sayin, S. (2014). A new algorithm for generating all nondominated solutions of multiobjective discrete optimization problems. European Journal of Operational Research, 232(3), 479–488. https://doi.org/10.1016/j.ejor.2013.08.001

[2] Tamby, S., & Vanderpooten, D. (2021). Enumeration of the nondominated set of multiobjective discrete optimization problems. INFORMS Journal on Computing, 33(1), 72–85. https://doi.org/10.1287/ijoc.2020.0953

[3] Mesquita-Cunha, M., Figueira, J. R., & Barbosa-Póvoa, A. P. (2022). New ϵ−constraint methods for multi-objective integer linear programming: A Pareto front representation approach. European Journal of Operational Research. https://doi.org/10.1016/j.ejor.2022.07.044
