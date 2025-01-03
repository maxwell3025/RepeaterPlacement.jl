"""
    CoordinateCostFunction(f, end_node_coordinate_matrix; kwargs...)
    CoordinateCostFunction(f, coords::Coordinates; kwargs...)
    CoordinateCostFunction(f, end_node_coordinate_matrix, kwargs, cache)

Cost function for repeater coordinates holding some constant parameters.

To call the cost function, either call `f` directly with a `Coordinates` object as
argument, or call the `CoordinateCostFunction` object as a function with a repeater
coordinate matrix as argument. It will automatically construct a `Coordinates` object
before calling `f`.
When calling this object directly, the `kwargs` are passed to `f` as keyword arguments.

The object has a `Dict{Symbol, Any}` stored under the `cache` field, which is meant to store
information about the latest cost-function evaluation.
Whenever the cost function is called, the `cache` is updated with the cost function value
under the key symbol `:cost_function_value`.
Whenever `f` is called, it is expected to return a `Dict{Symbol, Any}` as second value
which is then used to update the cache.
This allows storing and retrieving additional information about the result.
Conventionally, the `fill_cache` keyword argument can be passed with the value `true` to add
extra information to the cache that would not otherwise be included, but its implementation
depends on `f`.

While construction is possible by passing the `kwargs` and `cache` fields directly as
dictionaries, it is recommended to use the other two constructors, especially since
the cache should typically be empty upon initialization.

# Arguments
- `f`: the cost function to minimize. Should take a `Coordinates` object as only positional
    argument, and possibly a number of keyword arguments.
- `end_node_coordinate_matrix`: a matrix where each row represents the coordinates of an
    end node. If a `Coordinates` object is passed instead, the end-node coordinates are
    automatically extracted.
- `kwargs`: keyword arguments to pass to `f` when calling the `CoordinateCostFunction`
    object as a function.
"""
struct CoordinateCostFunction{F}
    f::F
    end_node_coordinate_matrix::Matrix{Float64}
    kwargs::Dict{Symbol, Any}
    cache::Dict{Symbol, Any}
end

function CoordinateCostFunction(f::F, end_node_coordinate_matrix; kwargs...) where F
    CoordinateCostFunction{F}(f, end_node_coordinate_matrix, Dict{Symbol, Any}(kwargs...),
        Dict{Symbol, Any}())
end

function CoordinateCostFunction(f, coords::Coordinates; kwargs...)
    CoordinateCostFunction(f, end_nodes(coords); kwargs...)
end

function (cost_fct::CoordinateCostFunction)(repeater_coordinate_matrix)
    coords = Coordinates(cost_fct.end_node_coordinate_matrix, repeater_coordinate_matrix)
    value, extra_cache = cost_fct.f(coords; cost_fct.cache...,  cost_fct.kwargs...)
    cost_fct.cache[:cost_function_value] = value
    for (key, v) in extra_cache
        cost_fct.cache[key] = v
    end
    value
end

Base.copy(x::CoordinateCostFunction) = CoordinateCostFunction(x.f,
    copy(x.end_node_coordinate_matrix), copy(x.kwargs), copy(x.cache))