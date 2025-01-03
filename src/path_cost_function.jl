"""
    pathwise_cost_matrix(end_node_coords, repeater_coords, cost_fct; radius=Inf,
        prune_fct=identity, max_num_paths=nothing, kwargs...)

Determine the cost associated with the best path for each pair of end nodes.

It is assumed that `cost_fct` is symmetric under reversion and that the available paths
from node i to node j are the same as those from node j to node i. Therefore, the resulting
cost matrix is symmetric.

# Arguments
- `coords`: a `Coordinates` object representing locations of end nodes and repeaters.
- `path_cost_fct`: a cost function that takes a `Path` object as argument.
- `radius`: the maximum distance between nodes for them to be connected.
    Having a small radius will allow for a faster calculation of the cost matrix as there
    will be fewer edges in the graph and hence fewer paths to consider. Moreover, it seems
    unlikely that the best paths will include very long links. However, if the radius is
    too small, the best possible path may not be considered.
- `prune_fct`: a function to prune the paths before determining the cost.
    This function should take a vector of paths and return a vector that contains a subset
    of the input paths. Only the returned paths are considered candidates for the best path.
    Having a selective prune function can result in a speedup by avoiding calculation of
    the cost function, but if it is too selective, it may prune out the best path.
- `max_num_paths`: the maximum number of paths to consider for each pair of end nodes.
    Will consider the shortest paths only. If `nothing`, all paths are considered.
    The pathfinding algorithm for finding the shortest paths is rather slow, so this will
    only grant a performance boost in pathfinding speed if the maximum number of paths is
    much smaller than the total number of paths. Additionally, it can function as a type of
    pruning, reducing the number of times the cost function needs to be evaluated.
- `kwargs...`: keyword arguments to be passed to `path_cost_fct`.
"""
function pathwise_cost_matrix(coords::Coordinates{T}, path_cost_fct;
        path_finding_function=path_finding_for_monotonic_costs, radius=Inf,
        linear_interpolation_num_nodes=10, verbose=false, kwargs...) where T
    g = build_graph(coords, radius)
    N = num_end_nodes(coords)
    cost_matrix = T(undef, N, N)
    cache_matrix = Matrix{Dict{Symbol, Any}}(undef, N, N)
    node_pairs = [(s, t) for s = 1:N for t = 1:(s-1)]
    for (s, t) in node_pairs
        g_internal = copy(g)
        r = copy(radius)
        guess_from_linear_interpolation =
            path_from_linear_interpolation(coords, s, t, linear_interpolation_num_nodes)
        node_sequence_guesses = [nodes(guess_from_linear_interpolation)]
        if :all_paths in keys(kwargs)
            path_guess = filter(kwargs[:all_paths]) do p
                (first(p.nodes) == s && last(p.nodes) == t) ||
                    (first(p.nodes) == t && last(p.nodes) == s)
            end
            path_guess = path_guess[1]
            push!(node_sequence_guesses, nodes(path_guess))
        end
        cost, cache = nothing, nothing
        while true
            path, cost, cache  = path_finding_function(path_cost_fct, g_internal, s, t;
                node_sequence_guesses=node_sequence_guesses, kwargs...)
            !isnothing(cost) && (cache[:path] = path; break)
            r += 10
            g_internal = build_graph(coords, r)
            verbose && println("No path found between nodes $i and $j for radius $r.
                Increasing radius by 10.")
        end
        cost_matrix[s, t] = cost_matrix[s, t] = cost
        cache_matrix[s, t] = cache_matrix[s, t] = cache
    end
    cost_matrix, cache_matrix
end

"""
    PathwiseCostFct

Abstract type for functions that assign cost to a `Coordinates` object based on a per-path
    cost between end nodes.

The per-path cost function (which takes a `Path` as its argument) can be retrieved using the
`path_cost_fct` function.
"""
abstract type PathwiseCostFct end

"""
    path_cost_fct(f::PathwiseCostFct)

Return the function defined on `Path` objects underlying a `PathwiseCostFct`.

A `PathWiseCostFct` is a function on `Coordinates` objects, but the value of the function
is based on an underlying per-path cost function.
This function returns that per-path cost function.
It can also be used on objects that wrap a `PathwiseCostFct`, such as a
`CoordinateCostFunction` or a `CoordinateSolution`.
"""
function path_cost_fct end

path_cost_fct(f::CoordinateCostFunction{F}) where F<:PathwiseCostFct = path_cost_fct(f.f)

"""
    largest_pathwise_cost(coords::Coordinates, path_cost_fct; radius=Inf,
        prune_fct=identity, max_num_paths=nothing, kwargs...)

Return cost associated to most expensive best path over all combination of end nodes.

This can be used as an optimization cost function for optimizing the worst cost over all
pairs of end nodes. This allows maximizing the guaranteed quality of service in the network.

# Arguments
- `coords`: a `Coordinates` object representing locations of end nodes and repeaters.
- `path_cost_fct`: a cost function that takes a `Path` object as argument.
- `resue_path`: if true, and if a `path` keyword argument is passed (which is typically
    passed from the cache of the previous solution), path finding is skipped and the cost
    is calculated assuming that the passed path is the best path corresponding to the worst
    pair of end nodes.
- `kwargs...`: keyword arguments to be passed to `path_cost_fct`.
"""
function largest_pathwise_cost(coords::Coordinates, path_cost_fct; reuse_path=false,
        kwargs...)
    if reuse_path && :path in keys(kwargs)
        path = Path(nodes(kwargs[:path]), build_graph(coords))
        cost, cache = path_cost_fct(path; kwargs...)
        cache[:path] = path
        return cost, cache
    end
    cost_mat, extra_cache_mat = pathwise_cost_matrix(coords, path_cost_fct; kwargs...)
    vectorized_cost = [cost_mat[i, j] for i = 1:num_end_nodes(coords) for j = 1:(i-1)]
    cost, index = findmax(vectorized_cost)
    vectorized_extra_cache = [extra_cache_mat[i, j] for i = 1:num_end_nodes(coords)
        for j = 1:(i-1)]
    extra_cache = vectorized_extra_cache[index]
    all_paths = [c[:path] for c in vectorized_extra_cache]
    extra_cache[:all_paths] = all_paths
    cost, extra_cache
end

"""
    LargestPathwiseCostFct(path_cost_fct)

Return a function-like object that determines the largest value of `path_cost_fct` over all
pairs of end nodes contained in a `Coordinates` object.

Any additional keyword arguments passed to the created function-like object are passed to
`path_cost_fct`.
The original path cost function is stored under the `path_cost_fct` field.

# Arguments
- `path_cost_fct`: a cost function that takes a `Path` object as argument.
"""
struct LargestPathwiseCostFct{F} <: PathwiseCostFct
    cost_fct::Function
    path_cost_fct::F
    function LargestPathwiseCostFct(path_cost_fct::F) where F
        f(coords::Coordinates; kwargs...) = largest_pathwise_cost(coords, path_cost_fct;
            kwargs...)
        new{F}(f, path_cost_fct)
    end
end

(f::LargestPathwiseCostFct)(coords::Coordinates; kwargs...) = f.cost_fct(coords; kwargs...)

path_cost_fct(f::LargestPathwiseCostFct) = f.path_cost_fct

function num_repeaters_on_paths(paths::Vector{Path{T}}) where T
    reps = Set()
    for path in paths
        path_reps = path.nodes[2:end-1]
        reps = union(reps, Set(path_reps))
    end
    length(reps)
end

"""
Mean cost, where mean is taken over the best paths for each pair of end nodes.

Note: currently untested and unused.
"""
function mean_pathwise_cost(coords::Coordinates, path_cost_fct; kwargs...)
    cost_mat, extra_cache_mat = pathwise_cost_matrix(coords, path_cost_fct; kwargs...)
    vectorized_cost = [cost_mat[i, j] for i = 1:num_end_nodes(coords) for j = 1:(i-1)]
    mean(vectorized_cost), Dict{Symbol, Any}()
    # currently no extra cache; unclear which cache to return when a mean is taken
end

"""
    MeanPathwiseCostFct(path_cost_fct)

Return a function-like object that determines the mean value of `path_cost_fct` over all
pairs of end nodes contained in a `Coordinates` object.
"""
struct MeanPathwiseCostFct{F} <: PathwiseCostFct
    cost_fct::Function
    path_cost_fct::F
    function MeanPathwiseCostFct(path_cost_fct::F) where F
        f(coords::Coordinates; kwargs...) = mean_pathwise_cost(coords, path_cost_fct;
            kwargs...)
        new{F}(f, path_cost_fct)
    end
end

(f::MeanPathwiseCostFct)(coords::Coordinates; kwargs...) = f.cost_fct(coords; kwargs...)

path_cost_fct(f::MeanPathwiseCostFct) = f.path_cost_fct