# RepeaterPlacement

If there are a number of end nodes in a Euclidean space that want to form a network together, and there is a limited number of repeater nodes that you can put anywhere in that space, where should you put them such that some network utility function is maximized?
This package has been developed to answer exactly that question!
With RepeaterPlacement, you can:
- Define the coordinates of end nodes and repeaters (using the `Coordinates` object);
- Define network utility functions (or cost functions) for coordinates;
- Use path-finding algorithms to find best paths through the network, allowing the definition of path-based network utility functions (where we assume that the path cost functions may not be isotonic, creatly complicating pathfinding as traditional methods like Dijkstra's algorithms fail);
- Find the best repeater coordinates using optimization algorithms based on gradients obtained through automatic differentiation;
- Visualize networks, paths and optimization results;
- Organize optimization data.

This package was written with quantum networks and quantum repeaters in mind, but there is nothing in this package that is specific to the quantum case.
Non-isotonic path cost functions are typical in quantum networks though, which has motivated the inclusion of pathfinding algorithms that work in that case.

For an example of how to use this package, it is recommended to check out the integration test included in the file `test/test_integration.jl`.


## Installation

Currently, RepeaterPlacement has not been added to the central Julia registry, as it is not yet clear to what extent this package will be used and developed further.
However, it can still be installed easily using the github url either by running
```julia
using Pkg; Pkg.add(url="https://github.com/GuusAvis/RepeaterPlacement.jl")
```
in Julia, or by using `]` in the REPL to go to Pkg mode and then typing
```julia
add https://github.com/GuusAvis/RepeaterPlacement.jl
```


## Citation

This package was developed in the process of writing the paper
Optimization of Quantum-Repeater Networks using Stochastic Automatic Differentiation
by Guus Avis and Stefan Krastanov.
Please refer to that paper when discussing or using this software package.

