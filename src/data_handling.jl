"""
    reformat_parameters(parameters::Dict;
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}())

Reformat a dictionary of parameters such that it is easier to save and interpret.

The start parameter and decay parameter of `kwargs[:schedule]` are extracted and added
to the output dictionary using the keys `:initial_temp` and `:decay_temp`.
This assumes the schedule is `Exp()` or similar from `ParameterSchedulers.jl`.

Each parameter that appears in `disallowed_kwargs` is not included in the output dictionary.
Furthermore, `kwarg_map` can be used to change the names of kwargs to make them more
readable or shorter.
"""
function reformat_parameters(parameters::Dict;
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}())
    if :schedule in keys(parameters)
        schedule = get(parameters, :schedule, nothing)
        parameters[:initial_temp] = schedule.start
        parameters[:decay_temp] = schedule.decay
    end
    d = Dict{Symbol, Union{Number, Base.UUID, PathEnumerationMethod}}()
    for (key, value) in parameters
        (key in disallowed_parameters || key == :schedule) && continue
        key = get(parameter_map, key, key)
        isnothing(value) && (value = 0)
        d[key] = value
    end
    d
end

"""
    get_savebase_name(sol::CoordinateSolution; digits=2,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}()), uuid=nothing)
    get_savebase_name(trace::CoordinateTrace; digits=2,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}()), uuid=nothing)

Determine the directory name to use for saving a solution.

The names and values (rounded according to the value of `digits`) of all the parameters in
the solution, processed by `reformat_parameters()`, are used.
`disallowed_parameters` and `parameter_map` are passed to this function call.
`uuid` is the unique identitifier used to save the solution.
If set to `nothing`, an attempt is made to obtain it from the solution.
If that fails, one is automatically generated.
"""
function get_savebase_name(sol::CoordinateSolution; digits=2,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(), uuid=nothing)
    number_of_reps = num_repeaters(Coordinates(sol))
    savebase = "n_reps_$(number_of_reps)"
    params = reformat_parameters(sol.X.kwargs,
        disallowed_parameters=disallowed_parameters, parameter_map=parameter_map)
    for (key, value) in params
        key == :uuid && continue
        if (value isa Number) && ! (value isa Bool)
            value = round(value, digits=digits)
        end
        savebase *= "_$(key)_$(value)"
    end
    if isnothing(uuid)
        if :uuid in keys(params)
            uuid = params[:uuid]
        else
            uuid = uuid4()
        end
    end
    savebase *= "-$(uuid)"
end

get_savebase_name(trace::CoordinateTrace; digits=2,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(), uuid=nothing) =
    get_savebase_name(trace[begin]; digits=digits,
        disallowed_parameters=disallowed_parameters, parameter_map=parameter_map, uuid=uuid)

"""
    save(sol::CoordinateSolution; digits=2, savebase=nothing,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}())
    save(trace::CoordinateTrace; digits=2, savebase=nothing,
        disallowed_parameters=[:tempearture, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}())

Save a solution or trace of solutions in a folder with name determined by `savebase`.

If `savebase` is `nothing`, one is generated automatically.
This includes a list of parameter names and values, rounded according to `digits`.
These are processed using `reformat_parameters()`.
`disallowed_parameters` and `parameter_map` are passed to this function call.

Solutions are saved using `JLD2`. A dictionary with the key `sol` is saved as `sol.jld2`
for a `CoordinateSolution`, a dictionary with the key `trace` is saved as `trace.jld2`
for a `CoordinateTrace`. The results can be easily loaded by calling `load` on the
directory name.
"""
function save(sol::CoordinateSolution; digits=2, savebase=nothing,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}())
    isnothing(savebase) && (savebase = 
        get_savebase_name(sol; digits=digits,
            disallowed_parameters=disallowed_parameters, parameter_map=parameter_map)
        )
    isdir(savebase) || mkdir(savebase)
    dict = Dict("sol" => sol)
    JLD2.save("$(savebase)/sol.jld2", dict)
end

# TODO move or remove function below
"""
    _get_processed_final_solution(trace::CoordinateTrace)

Get a version of `trace` that only has kwargs that occur in both first and last solutions.
"""
function _get_processed_final_solution(trace::CoordinateTrace)
    initial_sol = trace[begin]
    final_sol = copy(trace[end])
    # only include kwargs that appear in both initial and final solution
    filter!(x -> x in initial_sol.X.kwargs, final_sol.X.kwargs)
    final_sol
end

function save(trace::CoordinateTrace; digits=2, savebase=nothing,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(), kwargs...)
    dict = Dict("trace" => trace)
    isnothing(savebase) && (savebase = get_savebase_name(trace; digits=digits,
            disallowed_parameters=disallowed_parameters, parameter_map=parameter_map)
        )
    isdir(savebase) || mkdir(savebase)
    JLD2.save("$(savebase)/trace.jld2", dict)
end

"""
    load(directory::String)

Load a solution or trace of solutions from a directory created using `save()`.
"""
function load(directory::String)
    if isfile("$(directory)/trace.jld2")
        dict = JLD2.load("$(directory)/trace.jld2")
        return dict["trace"]
    end
    if isfile("$(directory)/sol.jld2")
        dict = JLD2.load("$(directory)/sol.jld2")
        return dict["sol"]
    end
    println("$directory contains no solution nor trace")
end

"""
    evaluate_complete(directory::String)

Load a saved `CoordinateTrace` and call `evaluate_complete()` on it. Save and return the results.

The evaluated results are saved in the same directory under the name "trace_evaluated.jld2".
"""
function evaluate_complete(directory::String; kwargs...)
    eval_file = "$(directory)/trace_evaluated.jld2"
    if isfile(eval_file)
        println("Skipping trace evaluation because this file already exists:\n
           $(eval_file)")
        return JLD2.load(eval_file)["trace"]
    end
    trace = load(directory)
    eval_trace = evaluate_complete(trace; kwargs...)
    dict = Dict("trace" => eval_trace)
    JLD2.save(eval_file, dict)
    eval_trace
end

"""
    has_cache(sol::CoordinateSolution)
    has_cache(sol::CoordinateTrace)

Check if a solution or trace has a nonempty cache.
"""
function has_cache(sol::CoordinateSolution)
    haskey(sol.X.kwargs, :fill_cache) && sol.X.kwargs[:fill_cache]
end

function has_cache(trace::CoordinateTrace)
    has_cache(trace[end]) & has_cache(trace[begin])
end

"""
    populate_cache(directory::String)

Populate the cache of saved data in case the cache is still empty. Then return data.
"""
function populate_cache(directory::String; kwargs...)
    data = load(directory)
    has_cache(data) && return data
    populate_cache!(data; kwargs...)
    save(data, savebase=directory)
    data
end

"""
    get_info_dict(sol::CoordinateSolution,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(), kwargs...)
    get_info_dict(sol::CoordinateTrace,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(), kwargs...)

Obtain a dictionary with information about a solution or trace of solutions.

Parameters are processed using `reformat_parameters()`, to which `disallowed_parameters`
and `parameter_map` are passed.

For solutions:
The parameter `:n_reps` is the total number of repeaters in a solution.
The parameter `:network_size` is the largest distance between any end nodes.
The parameter `:n_used_reps` is the number of repeaters that lie on a used path.

For traces:
The parameter `:n_iterations` is the total number of iterations of the optimization.
The parameter `:total_time` is the sum of the computational times taken for each iteration.
Parameters are furthermore added for the final solution in the trace.

Other parameters are obtained for the solution's kwargs and cache.
"""
function get_info_dict(sol::CoordinateSolution,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(); kwargs...)
    params = sol.X.kwargs
    s = evaluate_complete(sol; kwargs...)
    cache = filter(x -> last(x) isa Number, s.X.cache)
    info = merge(params, cache)
    coords = Coordinates(sol)
    info[:n_reps] = num_repeaters(coords)
    info[:network_size] = max_dist(end_nodes(coords))
    info[:n_used_reps] = num_repeaters_on_paths(s)
    reformat_parameters(info, disallowed_parameters=disallowed_parameters,
        parameter_map=parameter_map)
end

function get_info_dict(trace::CoordinateTrace,
        disallowed_parameters=[:temperature, :epoch, :fill_cache],
        parameter_map=Dict{Symbol, Symbol}(); kwargs...)
    num_iterations = length(trace)
    sol = _get_processed_final_solution(trace)
    info = get_info_dict(sol, disallowed_parameters, parameter_map; kwargs...)
    info[:total_time] = total_time(trace)
    info[:n_iterations] = num_iterations
    info
end

"""
    write_data_to_directory(directory)

Extract information from the data saved in `directory` and save it as `data.csv`.
"""
function write_data_to_directory(directory; kwargs...)
    data = load(directory)
    if isnothing(data)
        println("Could not write to this directory, because there is no data:\n $directory")
        return
    end
    info_dict = get_info_dict(data; kwargs...)
    df = DataFrame(info_dict)
    file = "$(directory)/data.csv"
    CSV.write(file, df)
end

"""
    save(fig::Plots.Plot, filename)

Save a figure as a PNG file.
"""
function save(fig::Plots.Plot, filename)
    savefig(fig, filename * ".png")
    savefig(fig, filename * ".svg")
end

"""
    save(fig::Plots.Anmiation, filename)

Save an anmiation as a GIF file.
"""
save(fig::Plots.Animation, filename) = gif(fig, filename * ".gif")

"""
    create_figures(directory::String)

Obtain data and evaluated data from a directory and use it to make and save figures.
"""
function create_figures(directory::String, extra_plot_functions=Function[]; kwargs...)
    data = populate_cache(directory; kwargs...)
    evaluated_data = evaluate_complete(directory; kwargs...)
    ps = create_figures(data, evaluated_data;
        extra_plot_functions=extra_plot_functions, kwargs...)
    for i in 1:length(ps)
        save(ps[i], "$(directory)/fig_$(i)")
    end
end

"""
    _process_data_single_directory(directory; figures=false, overwrite=false)

Process the data that was saved to a specific directory.

This function extracts and writes data using `write_data_to_directory()`.
Existing data is overwritten if and only if `overwrite` is set to `true`.
If `figures` is true, `create_figures` is also called.

If an EOF data is encountered when processing, this probably means that the data in the
directory is broken. This function will then delete that folder, so be careful!
"""
function _process_data_single_directory(directory; figures=false, overwrite=false,
        kwargs...)
    file_read_error = Union{EOFError, JLD2.InvalidDataException}
    try
        if overwrite || ! isfile(joinpath(directory, "data.csv"))
            write_data_to_directory(directory; kwargs...)
        end
        if figures && (overwrite || !isfile(joinpath(directory, "fig_1.png")))
            create_figures(directory; kwargs...)
        end
    catch e
        if typeof(e) isa file_read_error ||
                (hasfield(typeof(e), :ex) && e.ex isa file_read_error)
            println("Error processing directory: $directory")
            println("This probably means the files in the directory are defunct.")
            println("Will now delete the directory.")
            println("Error: $e")
            rm(directory; recursive=true)
        else
            println("Error processing directory: $directory")
            rethrow(e)
        end
    end
end

"""
    _process_data_many_directories(directory; figures=false, overwrite=false)

Process the data in all data directories contained by `directory`.

This function extracts and writes data using `write_data_to_directory()`.
Existing data is overwritten if and only if `overwrite` is set to `true`.
If `figures` is true, `create_figures` is also called.
"""
function _process_data_many_directories(directory; figures=false, overwrite=false,
        kwargs...)
    directories = filter(isdir, readdir(directory, sort=false, join=true))
    # TODO do we need the filter or can we simply catch errors?
    filter(x -> x[begin:5] == "n_reps", directories)
    for directory in directories
        _process_data_single_directory(directory; figures=figures, overwrite=overwrite,
            kwargs...)
    end
end

"""
    process_data(directory; batch=false, figures=false, overwrite=false)

Process data in a directory.

If `batch` is `false`, only the data that is contained directly in `directory` is processed.
If it is `true`, all data directories contained by `directory` are processed.

This function extracts and writes data using `write_data_to_directory()`.
Existing data is overwritten if and only if `overwrite` is set to `true`.
If `figures` is true, `create_figures` is also called.

If an EOF data is encountered when processing a data folder, this probably means that the
data in the directory is broken. This function will then delete that folder, so be careful!
"""
function process_data(directory; batch=false, figures=false, overwrite=false, kwargs...)
    if batch
        _process_data_many_directories(directory; figures=figures, overwrite=overwrite,
            kwargs...)
    else
        _process_data_single_directory(directory; figures=figures, overwrite=overwrite,
            kwargs...)
    end
end

function find_directory(uuid, directory=pwd())
    directories = filter(isdir, readdir(directory, sort=false, join=true))
    filter!(x->occursin(string(uuid), x), directories)
    directories[1]
end

function find_directory(sol::CoordinateSolution, directory=pwd())
    find_directory(sol.X.kwargs[:uuid], directory)
end

function find_directory(trace::CoordinateTrace, directory=pwd())
    find_directory(trace[end].X.kwargs[:uuid], directory)
end

"""
    resume_simulated_annealing(directory, num_epochs=nothing
        save::SaveType=DoNotSave, keep_old_file=true, verbose=false)

Resume a simulated-annealing optimization that has been saved to a directory.

After calling this function on a directory, the file `trace.jld2` stored therein will be
the result of the extended optimization, which is different from the original `trace.jld2`
in that directory.

# Parameters
- directory: the directory where the optimization was saved.
- num_epochs: the number of epochs to run the optimization for.
    If more than the original number of epochs, the optimization is extended.
- save: whether and how often to save the results of the optimization.
- keep_old_file: if true, the file storing the original trace is kept, but its name is
    changed from `trace.jld2` to `trace_<time_stamp>.jld2`.
- verbose: if true, print information about the optimization.
"""
function resume_simulated_annealing(directory, num_epochs=nothing;
        save::SaveType=DoNotSave, keep_old_file=true, verbose=false)
    trace = load(directory)
    save != DoNotSave && keep_old_file &&
        mv("$directory/trace.jld2", "$directory/trace_$(Dates.now()).jld2")
    resume_simulated_annealing!(trace; num_epochs=num_epochs, save=save, savebase=directory,
        verbose=verbose)
    trace
end