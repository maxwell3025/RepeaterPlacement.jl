"""
    StochasticAD.StochasticModel(X, p::Coordinates; kwargs...)

When a `Coordinates` is passed for `p`, a `CoordinateCostFunction` is automatically
constructed to create a `StochasticModel` that is also a `CoordinateSolution`.
Only the repeater coordinates are passed to the underlying `StochasticModel` to be
optimized.
"""
function StochasticAD.StochasticModel(X, p::Coordinates; kwargs...)
    cost_fct = CoordinateCostFunction(X, p; kwargs...)
    StochasticModel(cost_fct, p.repeater_coordinate_matrix)
end

"""
    StochasticAD.StochasticModel(X, p; kwargs...)

When keyword arguments are passed, these are disregarded.
This method was added for compatibility with optimization methods of `RepeaterPlacement`.
"""
StochasticAD.StochasticModel(X, p; kwargs...) = StochasticModel(X, p)

function Base.copy(x::StochasticModel)
    new_X = isbitstype(typeof(x.X)) ? x.X : copy(x.X)
    StochasticModel(new_X, copy(x.p))
end

"""
    CoordinateSolution(cost_function, coords::Coordinates; kwargs...)

A `StochasticAD.StochasticModel` with a `CoordinateCostFunction` as cost function.

This is the type used to parameterize solutions in the optimization of repeater coordinates.
A `Coordinates` object representing the solution can be obtained by calling `Coordinates` on
this object.
"""
CoordinateSolution{F, S} = StochasticModel{S, CoordinateCostFunction{F}}

Coordinates(x::CoordinateSolution) = Coordinates(x.X.end_node_coordinate_matrix, x.p)
num_end_nodes(x::CoordinateSolution) = num_end_nodes(Coordinates(x))
num_repeaters(x::CoordinateSolution) = num_repeaters(Coordinates(x))
num_nodes(x::CoordinateSolution) = num_nodes(Coordinates(x))
path_cost_fct(x::CoordinateSolution{F, S}) where {F<:PathwiseCostFct, S} =
    path_cost_fct(x.X)
num_repeaters_on_paths(x::CoordinateSolution{LargestPathwiseCostFct{F}, S}) where {F, S} =
    num_repeaters_on_paths(x.X.cache[:all_paths])

"""
    CoordinateTrace

Vector of `CoordinateSolution` objects. Used to represent the trace of optimizations.
"""
CoordinateTrace{F, S} = Vector{CoordinateSolution{F, S}}

total_time(trace::CoordinateTrace) = sum(sol.X.cache[:time] for sol in trace)

"""
    populate_cache!(sol::CoordinateSolution; fill_cache=true, kwargs...)
    populate_cache!(trace::CoordinateTrace; fill_cache=true, kwargs...)

Populate the cache of the `CoordinateCostFunction` stored as `sol.X`.

Any keyword arguments are passed on to the cost function.
Of particular note is the `fill_cache` keyword.
Cost functions can be set up to take the keyword `fill_cache`, such that they only evaluate
the cost as efficiently as possible when `fill_cache` is false but also calculate and store
additional information about the solution when `fill_cache` is true.
A typical workflow involves running an optimization using `fill_cache=false` (or
equivalently, not passing it to the optimizer), and then using `populate_cache!` with
`fill_cache=true` when inspecting results.

Values that are already in the cache may be overwritten by this function.
"""
function populate_cache!(sol::CoordinateSolution; fill_cache=true, kwargs...)
    sol.X.kwargs[:fill_cache] = fill_cache
    for (key, value) in kwargs
        sol.X.kwargs[key] = value
    end
    sol.X(sol.p)
end

function populate_cache!(trace::CoordinateTrace; fill_cache=true, kwargs...)
    for sol in trace
        populate_cache!(sol; fill_cache=fill_cache, kwargs...)
    end
end

"""
    populate_cache(sol; fill_cache=true, kwargs...)

Non-mutating version of `populate_cache!`.
"""
function populate_cache(sol; fill_cache=true, kwargs...)
    s = copy(sol)
    populate_cache!(s; fill_cache=fill_cache, kwargs...)
    s
end

"""
    evaluate_complete(sol::CoordinateSolution; kwargs...)
    evaluate_complete(trace::CoordinateTrace; kwargs...)

Return a copy of `sol` with the cache filled using zero temperature and whatever additional
keyword arguments are passed.

Sometimes the evaluation of a cost function is noisy, with noise levels quantified by
the `temperature` keyword argument (e.g., by tuning the number of samples taken in a
Monte-Carlo simulation). This keyword argument is used by the
`RepeaterPlacement.simulated_annealing()` optimization function.

This function creates a copy of `sol` and calls `populate_cache!` on it with the keyword
arguments `populate_cache=true` and `tempearture=0`, resulting in the most complete and
noiseless results that can be obtained. Often useful when inspecting optimization results,
e.g. before calling `plot_cost_trace()`. See also the documentation of `populate_cache!()`.

Note that `evaluate_complete()` may be slow compared to direct evaluation of the cost
function without populating the cache and setting the temperature to zero.

If called on a `CoordinateTrace`, each solution in the trace will be evaluated completely.
"""
function evaluate_complete(sol::CoordinateSolution; kwargs...)
    s = copy(sol)
    populate_cache!(s; fill_cache=true, temperature=0., kwargs...)
    s
end

function evaluate_complete(trace::CoordinateTrace, args...; kwargs...)
    [evaluate_complete(sol, args...; kwargs...) for sol in trace]
end

"""
    gradient_descent(f, initial_sol, num_iterations, learning_rate=0.1; kwargs...)

Perform a gradient descent to (local) minimum of function `f`.

Return trace of solutions as a vector. The first element of the trace is the
initial solution `initial_sol`, the final element of the trace is the final solution.

This function was written with `CoordinateCostFunction` in mind for the input `f`,
but it should work with other function-like objects as well.
If `f` has a field `:kwargs` that holds a dictionary, as `CoordinateCostFunction` does,
any additional keyword arguments (`kwargs...`) are written here.
For `CoordinateCostFunction`, this guarantees that these keyword arguments are passed
to the cost function in each evaluation.
The parameters `num_iterations` and `learning_rate` are also stored in the `kwargs`
dictionary of the cost function.

# Arguments

- `f`: the function to minimize.
- `initial_sol`: the initial solution to start the gradient descent from.
- `num_iterations`: the number of iterations to perform.
- `learning_rate`: the learning rate used for the Adam optimizer.
"""
function gradient_descent(f, initial_sol, num_iterations, learning_rate=0.1; kwargs...)
    model = StochasticModel(f, copy(initial_sol); num_iterations=num_iterations,
        learning_rate=learning_rate, kwargs...)
    cost_fct_has_cache = hasfield(typeof(model.X), :cache)
    optimiser = Adam(learning_rate)
    s = Optimisers.setup(optimiser, model)
    initial_model = copy(model)
    if cost_fct_has_cache
        populate_cache!(initial_model; fill_cache=false, kwargs...)
        initial_model.X.cache[:time] = 0.
    end
    trace = typeof(model)[initial_model]
    for _ in 1:num_iterations
        time = update_and_time!(s, model)
        new_model = copy(model)
        cost_fct_has_cache && (new_model.X.cache[:time] = time)
        push!(trace, new_model)
    end
    trace
end

"""
How often to save the trace during optimization.
"""
@enum SaveType DoNotSave=1 SaveAtTheEnd=2 SaveEveryEpoch=3 SaveEveryIteration=4

"""
    simulated_annealing(f, initial_sol, num_epochs, epoch_size;
        schedule=Exp(start=10., decay=0.8), verbose=false, save::SaveType=DoNotSave,
        savebase=nothing, kwargs...)

Perform simulated annealing on the learning rate of the Adam optimizer and return trace.

Function `f` is optimized using stochastic-gradient optimizer Adam, where a "temperature"
is used to set Adam's learning rate which is gradually decreased over the optimization.
The optimization is done through `num_epochs` epochs, each consisting of `epoch_size` steps.
Each epoch has a fixed temperature, which is determined by the `schedule` function.
The function `f` is also passed the temperature as a keyword argument, such that further
temperature-dependent behavior can be implemented in `f`.

This function was written with `CoordinateCostFunction` in mind for the input `f`,
but it should work with other function-like objects as well.
If `f` has a field `:kwargs` that holds a dictionary, as `CoordinateCostFunction` does,
any additional keyword arguments (`kwargs...`) are written here.
For `CoordinateCostFunction`, this guarantees that these keyword arguments are passed
to the cost function in each evaluation.
The parameters `num_epochs`, `epoch_size` and `schedule` are also stored in the `kwargs`
dictionary of the cost function, along with the UUID used to store the results in case
data is saved.
Additionally, during the optimization, the current epoch number and temperature are
stored in the `kwargs` dictionary.
Automated saving is currently only supported for `CoordinateCostFunction`, to make it work
for other types of cost functions the appropriate methods of `RepeaterPlacement.save`
should be implemented.

# Arguments
- `f`: the function to minimize.
- `initial_sol`: the initial solution to start the optimization from.
- `num_epochs`: the total number of epochs to complete in simulated annealing.
- `epoch_size`: the number of iterations per epoch.
- `schedule`: the schedule function (probably a `ParameterSchedulers.AbstractSchedule`)
    to use in the optimization.
- `verbose`: whether to print progress information.
- `save::SaveType`: indicates if and how often to save the trace during optimization.
- `savebase`: the base name to use for saving the trace.
    If `nothing`, a new UUID is generated and used as the base name.
    Ignored if `save` is `DoNotSave`.
"""
function simulated_annealing(f, initial_sol, num_epochs, epoch_size;
        schedule=Exp(start=10., decay=0.8), verbose=false, save::SaveType=DoNotSave,
        savebase=nothing, kwargs...)
    model = StochasticModel(f, copy(initial_sol); temperature=first(schedule), epoch=1,
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=schedule, kwargs...)
    cost_fct_has_cache = hasfield(typeof(model.X), :cache)
    if cost_fct_has_cache
        populate_cache!(model; fill_cache=false, kwargs...)
        model.X.cache[:time] = 0.
    end
    trace = typeof(model)[model]
    resume_simulated_annealing!(trace, num_epochs=num_epochs, epoch_size=epoch_size,
        schedule=schedule, verbose=verbose, save=save, savebase=savebase)
end

"""
    resume_simulated_annealing!(trace; num_epochs=nothing, epoch_size=nothing,
        schedule=nothing, verbose=false, save::SaveType=DoNotSave, savebase=nothing)

Resume simulated annealing optimization from a trace.

This function can be used to continue an optimization that was interrupted but saved,
or to increase the number of iterations of an optimization that has already been completed.

See also the documtentation of `simulated_annealing()`.

# Arguments
- trace: a `Vector{StochasticAD.StochasticModel}` of solutions so far that should be
    extended.
- num_epochs: the total number of epochs to complete in simulated annealing.
    Should be at least as larger as the number of epochs already completed.
    If `nothing`, and if the solutions in the trace carry their own value for the number
    of epochs, the value from the trace is used.
    Note that the value carried by the trace may not be equal to the number of epochs that
    were actually completed to create the trace, e.g., when the optimization crashes
    before finishing but the results obtained so far were saved.
    In that case, this function can be used to continue that optimization.
    If a value is provided and the trace carries a smaller value, the value in the trace
    is updated.
    This can be used to extend the number of epochs in the optimization beyond
    the originally set value.
- epoch_size: the number of iterations per epoch. Should be equal to the epoch size in the
    original optimization, although this is only enforced if the solutions in the trace
    carry their own value for the epoch size.
    If `nothing`, the value is automatically extracted from the trace if possible.
- schedule: the schedule function (probably a `ParameterSchedulers.AbstractSchedule`)
    to use in the optimization. Should be equal to the schedule used in the original
    optimization, although this is only enforced if the solutions in the trace carry their
    own schedule.
    If `nothing`, the schedule is extracted from the trace if possible.
- verbose: whether to print progress information.
- `save::SaveType`: indicates if and how often to save the trace during optimization.
- `savebase`: the base name to use for saving the trace.
    If `nothing`, a new UUID is generated and used as the base name.
    Ignored if `save` is `DoNotSave`.
"""
function resume_simulated_annealing!(trace; num_epochs=nothing, epoch_size=nothing,
        schedule=nothing, verbose=false, save::SaveType=DoNotSave, savebase=nothing)

    first_iteration_to_do = length(trace)
    active = isone(first_iteration_to_do) ? true : false
    !active && (current_iteration = 0)

    model = copy(trace[end])
    cost_fct_has_cache = hasfield(typeof(model.X), :cache)
    cost_fct_has_kwargs = hasfield(typeof(model.X), :kwargs)

    if cost_fct_has_kwargs && :schedule in keys(model.X.kwargs)
        if isnothing(schedule)
            schedule = model.X.kwargs[:schedule]
        else
            schedule == model.X.kwargs[:schedule] ||
                throw(
                    ArgumentError(
                        "passed schedule does not match schedule stored in trace"
                    )
                )
        end
    end
    isnothing(schedule) && throw(ArgumentError("schedule not provided"))

    if cost_fct_has_kwargs && :epoch_size in keys(model.X.kwargs)
        if isnothing(epoch_size)
            epoch_size = model.X.kwargs[:epoch_size]
        else
            epoch_size == model.X.kwargs[:epoch_size] ||
                throw(
                    ArgumentError(
                        "passed epoch_size does not match epoch_size stored in trace"
                    )
                )
        end
    end
    isnothing(epoch_size) && throw(ArgumentError("epoch_size not provided"))

    if cost_fct_has_kwargs && :num_epochs in keys(model.X.kwargs)
        if isnothing(num_epochs)
            num_epochs = model.X.kwargs[:num_epochs]
        else
            for sol in trace
                sol.X.kwargs[:num_epochs] = num_epochs
                model.X.kwargs[:num_epochs] = num_epochs
            end
        end
    end
    isnothing(num_epochs) && throw(ArgumentError("num_epochs not provided"))

    if save != DoNotSave && isnothing(savebase)
        if cost_fct_has_kwargs && :uuid in keys(model.X.kwargs)
            uuid = model.X.kwargs[:uuid]
        else
            uuid = uuid4()
            if cost_fct_has_kwargs
                for sol in trace
                    sol.X.kwargs[:uuid] = uuid
                end
                model.X.kwargs[:uuid] = uuid
            end
        end
        savebase = get_savebase_name(trace, uuid=uuid)
    end

    optimiser = Adam()
    opt_st = Optimisers.setup(optimiser, model)
    for e in 1:num_epochs
        temperature = schedule(e)
        verbose && active && println("starting epoch $e with temperature $temperature")
        if cost_fct_has_kwargs
            model.X.kwargs[:temperature] = temperature
            model.X.kwargs[:epoch] = e
        end
        Optimisers.adjust!(opt_st, temperature)  # use temperature as learning rate
        for _ in 1:epoch_size
            if !active
                current_iteration += 1
                if current_iteration == first_iteration_to_do
                    active = true
                    verbose && println("resuming at epoch $e with temperature $temperature \
                        after skipping the first $first_iteration_to_do iterations")
                end
            end
            active || continue
            time = update_and_time!(opt_st, model)
            new_model = copy(model)
            cost_fct_has_cache && (new_model.X.cache[:time] = time)
            push!(trace, new_model)
            save == SaveEveryIteration && RepeaterPlacement.save(trace; savebase=savebase)
        end
        save == SaveEveryEpoch && RepeaterPlacement.save(trace; savebase=savebase)
    end
    save == SaveAtTheEnd && RepeaterPlacement.save(trace; savebase=savebase)
    trace
end

"""
    update_and_time!(opt_st, model)

Update the model using the optimizer `opt_st` and return the time it took to update.

This function is used in the update step of the inner loops of optimization functions.
When a `LargestPathwiseCostFct` is being optimized, this function takes a shortcut
by first choosing a path to use, and then only determine the derivative for that fixed path.
"""
function update_and_time!(opt_st, model)
    @elapsed Optimisers.update!(opt_st, model, stochastic_gradient(model))
end

function update_and_time!(opt_st, model::CoordinateSolution{LargestPathwiseCostFct{S}, T}
        where {S, T})
    t1 = @elapsed model.X(model.p)
    model.X.kwargs[:reuse_path] = true
    t2 = @elapsed Optimisers.update!(opt_st, model, stochastic_gradient(model))
    model.X.kwargs[:reuse_path] = false
    t1 + t2
end