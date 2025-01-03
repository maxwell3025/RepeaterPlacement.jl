"""
This integration test also serves as an example of how to use this package.
It can be run as an interactive notebook using, e.g., VS Code with the Julia extension.
Make sure though to first run the "using" statements at the top of `runtests.jl`.
"""

## initialize
path_cost(p::Path; kwargs...) =
    sum(l^2 for l in p.lengths), Dict{Symbol, Any}(:num_reps => length(p) - 1)
cost_fct = LargestPathwiseCostFct(path_cost)
initial_solution = Coordinates()
initial_solution = add_end_node(initial_solution, [0., 0.])
initial_solution = add_end_node(initial_solution, [100., 0.])
initial_solution = add_repeater(initial_solution, [80., 0.])
p1 = plot_graph(initial_solution)  # show nodes and their possible connections as a graph
p2 = plot_node_locations(initial_solution)  # show nodes geometrically
@test p2 isa Plots.Plot

## solve using gradient descent
num_iterations = 50
learning_rate = 1.
trace = RepeaterPlacement.gradient_descent(cost_fct, initial_solution, num_iterations,
    learning_rate)
@test length(trace) == num_iterations + 1
sol = copy(trace[end])
@test 48 < first(repeater(Coordinates(sol), 1)) < 52
@test total_time(trace) > 0
for s in trace
    @test :time in keys(s.X.cache)
    @test :cost_function_value in keys(s.X.cache)
    @test Set(s.X.cache[:path].nodes) == Set([1, 2, 3])
    @test s.X.cache[:all_paths][1] == s.X.cache[:path]
end
p1 = plot_cost_trace(trace)
@test p1 isa Plots.Plot
p2 = plot_cache_value(trace, :num_reps, show_temperature=false)
@test p2 isa Plots.Plot
p3 = plot_node_locations(sol)  # show final solution
@test p3 isa Plots.Plot
p4 = plot_node_locations(trace, show_paths=true)  # show animation of the optimization
@test p4 isa Plots.Animation
ps = RepeaterPlacement.create_figures(trace)  # make bunch of figures in one go
@test length(ps) > 1

## solve using simulated annealing
num_epochs = 5
epoch_size = 5 
schedule = Exp(start=2., decay=0.8)
trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
    epoch_size; schedule=schedule, save=RepeaterPlacement.SaveEveryIteration,
    method=WavefrontSearch())
@test length(trace) == num_epochs * epoch_size + 1
sol = copy(trace[end])
@test 48 < first(repeater(Coordinates(sol), 1)) < 52
@test total_time(trace) > 0
for s in trace
    @test :time in keys(s.X.cache)
    @test :cost_function_value in keys(s.X.cache)
    @test Set(s.X.cache[:path].nodes) == Set([1, 2, 3])
    @test s.X.cache[:all_paths][1] == s.X.cache[:path]
    @test :temperature in keys(s.X.kwargs)
    @test 0 < s.X.kwargs[:temperature] ≤ 2.
    @test :epoch in keys(s.X.kwargs)
    @test 1 ≤ s.X.kwargs[:epoch] ≤ 5
end
p1 = plot_cost_trace(trace, show_temperature=true)
@test p1 isa Plots.Plot
p2 = plot_cache_value(trace, :num_reps, show_temperature=false)
@test p2 isa Plots.Plot
p3 = plot_node_locations(sol)  # show final solution
@test p3 isa Plots.Plot
p4 = plot_node_locations(trace)  # show animation of the optimization
@test p4 isa Plots.Animation
ps = RepeaterPlacement.create_figures(trace)  # make bunch of figures in one go
@test length(ps) > 1

## test saving, loading and data processing
save_dir = RepeaterPlacement.find_directory(sol)
trace_loaded = RepeaterPlacement.load(save_dir)
@test length(trace_loaded) == length(trace)
@test total_time(trace_loaded) == total_time(trace)
@test :schedule in keys(trace_loaded[begin].X.kwargs)  # there was a bug removing :schedule 
@test :uuid in keys(trace_loaded[begin].X.kwargs)  # there was a bug removing :uuid
sol_loaded = trace_loaded[end]
@test sol_loaded.p == sol.p
@test Coordinates(sol_loaded).end_node_coordinate_matrix ==
    Coordinates(sol).end_node_coordinate_matrix
RepeaterPlacement.process_data(save_dir, figures=true)
@test isfile(save_dir * "/data.csv")
@test isfile(save_dir * "/trace.jld2")
@test isfile(save_dir * "/trace_evaluated.jld2")
@test isfile(save_dir * "/fig_1.png")
@test isfile(save_dir * "/fig_2.gif")
@test isfile(save_dir * "/fig_3.png")

# test resuming an optimization to extend results
trace_files = filter(f -> startswith(f, "trace"), readdir(save_dir))
@test length(trace_files) == 2  # trace.jld2 and trace_evaluated.jld2
new_trace = RepeaterPlacement.resume_simulated_annealing(save_dir, num_epochs + 1,
    save=RepeaterPlacement.SaveAtTheEnd, keep_old_file=true)
new_trace_loaded = RepeaterPlacement.load(save_dir)
@test length(new_trace) == length(trace) + epoch_size
@test length(new_trace_loaded) == length(new_trace)
@test isfile(save_dir * "/trace.jld2")
trace_files = filter(f -> startswith(f, "trace"), readdir(save_dir))
@test length(trace_files) == 3  # the old trace.jld2 has been archived with a date
newer_trace = RepeaterPlacement.resume_simulated_annealing(save_dir, num_epochs + 2,
    save=RepeaterPlacement.SaveAtTheEnd, keep_old_file=false)
@test length(newer_trace) == length(new_trace) + epoch_size
trace_files = filter(f -> startswith(f, "trace"), readdir(save_dir))
@test length(trace_files) == 3  # the old trace.jld2 has been replaced
newer_trace_loaded = RepeaterPlacement.load(save_dir)
@test length(newer_trace_loaded) == length(newer_trace)

rm(save_dir, recursive=true)