# RAMODO
Reformulation Approaches for Multi-objective Discrete Optimization

This repo is created to include the supporting material for the TUBİTAK project 122M422 - "Reformulation Approaches for Multi-objective Discrete Optimization Problems" and will not be maintained in the future. To solve multiobjective optimization problems, use [MultiObjectiveAlgorithms.jl](https://github.com/jump-dev/MultiObjectiveAlgorithms.jl).

You can load problems and solve them using the algorithms provided.

First, run the main script `MODO.jl`. After that, load the problem using the data structures provided. Code works with the formats given in `txt` folder for each problem set.

```
include("MODO.jl")
```
$$
\begin{align*}
\max & \sum_{i=1}^p v_{i}^k x_i, \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n w_i x_i \le W \\
& x_i \in \{0, 1\}
\end{align*}
$$

## Knapsack Problem
```
kp = MOKnapsackProblem("Kirlik & Sayın/KP/txt/KP_p-3_n-10_ins-1.txt")
```
$$
\begin{align*}
\max & \sum_{i=1}^p c_{ij}^k x_{ij}, & \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n x_{ij} = 1, & \quad \forall j \in 1,...,n \\
& \sum_{j=1}^n x_{ij} = 1, & \quad \forall i \in 1,...,n \\
& x_{ij} \in \{0, 1\}
\end{align*}
$$

## Assignment Problem
```
ap = MOAssignmentProblem("Kirlik & Sayın/AP/txt/AP_p-3_n-5_ins-1.txt")
```
$$
\begin{align*}
\max & \sum_{i=1}^p c_{i}^k x_i, & \quad k=1,...,p \\
\text{s.to:} & \sum_{i=1}^n a_{ji} x_i \le b_j, & \quad \forall j \in 1,...,m \\
& x_i \in \{0, 1\}
\end{align*}
$$

## Integer Linear Problem
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

1) Overall Nondominated Vector Generation Ratio (ONVGR)
$$\text{ONVGR} = {|\mathcal{Y}_R| \over |\mathcal{Y}_N|}$$
`ongvr(Y_R, Y_N)`
2) Additive Epsilon Indicator ($\varepsilon_+$)
$$\varepsilon_+ = \max_{x \in \mathcal{Y}_N} \min_{y \in \mathcal{Y}_R} \max_{i = 1,...,p} {y_i - x_i \over r_i}$$
where $r_i$ is the range of objective $i$ in $\mathcal{Y}_N$
`ε₊(Y_R, Y_N)`
3) Coverage Error
$$\varepsilon_+ = \max_{x \in \mathcal{Y}_N} \min_{y \in \mathcal{Y}_R} \max_{i = 1,...,p} {|y_i - x_i| \over r_i}$$
where $r_i$ is the range of objective $i$ in $\mathcal{Y}_N$
`coverage_error(Y_R, Y_N)``
4) Range Ratio
$$\max_{i \in \{1, \dots, p\}}{\max_{y_i \in \mathcal{Y}_R} y_i - \min_{y_i \in \mathcal{Y}_R} y_i \over \max_{y_i \in \mathcal{Y}_N} y_i - \min_{y_i \in \mathcal{Y}_N} y_i}$$
`range_ratio(Y_R, Y_N)``
5) Hypervolume
See the [link](https://lopez-ibanez.eu/hypervolume) for computing $hv(\mathcal{Y}_R)$
`hypervolume(Y_R, reference)`
In order to use julia function to compute hypervolume, you first need to build a shared library (extensions are `.so` for Linux, `.dll` for Windows and `.dylib`for Mac) using the code provided with link.

6) Hypervolume Ratio
$${hv(\mathcal{Y}_R) \over hv(\mathcal{Y}_N)}$$
