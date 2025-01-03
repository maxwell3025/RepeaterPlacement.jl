@testset "CoordinateSolution and CoordinateTrace" begin
    # make a CoordinateSolution and test its coordinate-based properties
    cost_fct = (c; kwargs...) -> (num_nodes(c), Dict{Symbol, Any}(kwargs...))
    coordinates = Coordinates()
    coordinates = add_end_node(coordinates, [0., 0.])
    coordinates = add_end_node(coordinates, [10., 0.])
    coordinates = add_repeater(coordinates, [5., 0.])
    sol = StochasticAD.StochasticModel(cost_fct, coordinates, test="test")
    @test typeof(sol) == CoordinateSolution{typeof(cost_fct), Matrix{Float64}}
    @test coordinates == Coordinates(sol)
    @test num_end_nodes(sol) == num_end_nodes(coordinates)
    @test num_repeaters(sol) == num_repeaters(coordinates)
    @test num_nodes(sol) == num_nodes(coordinates)

    # evaluate the function (populates the cache)
    @test sol.X(sol.p) == 3
    @test sol.X.kwargs == Dict(:test => "test")
    @test sol.X.cache == Dict(:cost_function_value => 3, :test => "test")

    # test copy
    sol2 = copy(sol)
    @test Coordinates(sol2).end_node_coordinate_matrix ==
        coordinates.end_node_coordinate_matrix
    @test Coordinates(sol2).repeater_coordinate_matrix ==
        coordinates.repeater_coordinate_matrix
    @test sol.X.cache == sol2.X.cache
    @test sol.X.kwargs == sol2.X.kwargs
    @test sol.p == sol2.p

    # populate cache
    # first test if it overwrites existing values in kwargs and cache
    populate_cache!(sol, test="test2")
    @test sol.X.kwargs == Dict(:test => "test2", :fill_cache => true)
    @test sol.X.cache ==
        Dict(:cost_function_value => 3, :test => "test2", :fill_cache => true)
    # also test that it works on a fresh solution
    sol = StochasticAD.StochasticModel(cost_fct, coordinates, test="test")
    sol2 = populate_cache(sol, test="test2")
    populate_cache!(sol, test="test2")
    for s in [sol, sol2]
        @test s.X.kwargs == Dict(:test => "test2", :fill_cache => true)
        @test s.X.cache ==
            Dict(:cost_function_value => 3, :test => "test2", :fill_cache => true)
    end
    @test !(sol === sol2)
    
    # evaluate completely
    sol = StochasticAD.StochasticModel(cost_fct, coordinates, test="test")
    sol = evaluate_complete(sol, test="test3")
    @test sol.X.kwargs == Dict(:test => "test3", :fill_cache => true,
        :temperature => 0.)
    @test sol.X.cache ==
        Dict(:cost_function_value => 3, :test => "test3", :fill_cache => true,
            :temperature => 0.)

    # test CoordinateTrace
    sol = StochasticAD.StochasticModel(cost_fct, coordinates, test="test")
    tr = [sol, copy(sol)]
    @test typeof(tr) == CoordinateTrace{typeof(cost_fct), Matrix{Float64}}
    populate_cache!(tr, test="test2")
    for s in tr
        @test sol.X.kwargs[:test] == "test2"
        @test sol.X.cache[:test] == "test2"
        @test sol.X.cache[:cost_function_value] == 3
    end
    tr[1].X.cache[:time] = 1
    tr[2].X.cache[:time] = 10
    @test total_time(tr) == 11
end

@testset "gradient descent" begin
    cost_fct(x) = (x[begin] - 5.)^2 + 10.
    initial_solution = [10.]
    num_iterations = 20
    learning_rate = 1.
    trace = RepeaterPlacement.gradient_descent(cost_fct, initial_solution, num_iterations,
        learning_rate)
    @test length(trace) == num_iterations + 1
    for sol in trace
        @test sol.X == cost_fct
        @test 0. < sol.p[begin] < 11
    end
    @test 3. < trace[end].p[begin] < 7.
    @test !(trace[begin].p[begin] ≈ trace[end].p[begin])
end

# extend equality test between stochastic models for easier testing
import Base.==
==(m1::StochasticModel, m2::StochasticModel) =
    getfield.(Ref(m1), fieldnames(StochasticModel)) ==
        getfield.(Ref(m2), fieldnames(StochasticModel))

@testset "simulated annealing" begin
    cost_fct(x) = (x[begin] - 5.)^2 + 10.
    initial_solution = [10.]
    num_epochs = 6
    epoch_size = 5
    schedule = Exp(start=10., decay=0.8)
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule)
    @test length(trace) == num_epochs * epoch_size + 1
    for sol in trace
        @test sol.X == cost_fct
        @test 0. < sol.p[begin] < 11
    end
    @test 3. < trace[end].p[begin] < 7.
    @test !(trace[begin].p[begin] ≈ trace[end].p[begin])

    # test resuming an optimization
    trace_before_resume = deepcopy(trace)
    trace_after_empty_resume = deepcopy(trace_before_resume)
    # extending with same parameters should not do anything
    RepeaterPlacement.resume_simulated_annealing!(trace_after_empty_resume;
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=schedule)
    @test length(trace_before_resume) == length(trace_after_empty_resume)
    @test trace_before_resume[end] == trace_after_empty_resume[end]
    @test all(trace_before_resume .== trace_after_empty_resume)
    trace_after_extension = deepcopy(trace_before_resume)
    # this should add one extra epoch to the trace
    RepeaterPlacement.resume_simulated_annealing!(trace_after_extension;
        num_epochs=num_epochs + 1, epoch_size=epoch_size, schedule=schedule)
    @test length(trace_after_extension) == length(trace_before_resume) + epoch_size
    # if any of the keyword arguments are not supplied, an error should be thrown
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume);
        num_epochs=num_epochs, epoch_size=nothing, schedule=schedule)
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume);
        num_epochs=nothing, epoch_size=epoch_size, schedule=schedule)
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume);
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=nothing)
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume))
    # simulate an optimization stopping in the middle of an epoch
    trace_unfinished = deepcopy(trace_before_resume[1:end - floor(Int, epoch_size / 2)])
    trace_unfinished_resumed = deepcopy(trace_unfinished)
    RepeaterPlacement.resume_simulated_annealing!(trace_unfinished_resumed;
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=schedule)
    @test length(trace_unfinished) < length(trace_unfinished_resumed)
    @test length(trace_unfinished_resumed) == length(trace_before_resume)
    @test trace_unfinished[end].p[begin] ≠ trace_unfinished_resumed[end].p[begin]
end

struct CostFunctionWithCacheAndKwargs
    cache::Dict{Symbol, Any}
    kwargs::Dict{Symbol, Any}
end
Base.copy(c::CostFunctionWithCacheAndKwargs) =
    CostFunctionWithCacheAndKwargs(copy(c.cache), copy(c.kwargs))
function (c::CostFunctionWithCacheAndKwargs)(x)
    c.cache[:cache_from_evaluation] = true
    (x[begin] - 5.)^2 + 10.
end
function RepeaterPlacement.populate_cache!(
        m::StochasticModel{S, CostFunctionWithCacheAndKwargs}; kwargs...) where S
    m.X.cache[:cache_populated] = true
    for (key, value) in kwargs
        m.X.kwargs[key] = value
    end
    m.X(m.p)
end
==(c1::CostFunctionWithCacheAndKwargs, c2::CostFunctionWithCacheAndKwargs) =
    c1.kwargs == c2.kwargs && c1.cache == c2.cache
function StochasticAD.StochasticModel(X::CostFunctionWithCacheAndKwargs, p::S;
        kwargs...) where S<:AbstractArray
    merge!(X.kwargs, kwargs)
    @invoke StochasticModel(X::Any, p::S)
end
saved = 0
function RepeaterPlacement.save(trace::Vector{StochasticModel{S,
        CostFunctionWithCacheAndKwargs}}; savebase) where S
    @test savebase == "test"
    global saved
    saved += 1
end

@testset "simulated annealing with cache and kwargs" begin
    cost_fct = CostFunctionWithCacheAndKwargs(Dict{Symbol, Any}(), Dict{Symbol, Any}())
    @test hasfield(typeof(cost_fct), :cache)
    @test hasfield(typeof(cost_fct), :kwargs)
    initial_solution = [10.]
    num_epochs = 6
    epoch_size = 5
    schedule = Exp(start=10., decay=0.8)
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule, test_kwarg="test")
    @test length(trace) == num_epochs * epoch_size + 1
    for (i, sol) in enumerate(trace)
        @test typeof(sol.X) == CostFunctionWithCacheAndKwargs
        @test 0. < sol.p[begin] < 11
        @test sol.X.cache[:cache_populated] == true
        @test sol.X.cache[:time] ≥ 0.
        @test sol.X.kwargs[:num_epochs] == num_epochs
        @test sol.X.kwargs[:epoch_size] == epoch_size
        @test sol.X.kwargs[:schedule] == schedule
        @test sol.X.kwargs[:temperature] ≥ 0.
        @test 1 ≤ sol.X.kwargs[:epoch] ≤ num_epochs
        @test sol.X.kwargs[:test_kwarg] == "test"
        @test sol.X.cache[:cache_from_evaluation] == true
    end
    @test 3. < trace[end].p[begin] < 7.
    @test !(trace[begin].p[begin] ≈ trace[end].p[begin])

    trace_before_resume = deepcopy(trace)
    trace_after_empty_resume = deepcopy(trace_before_resume)
    RepeaterPlacement.resume_simulated_annealing!(trace_after_empty_resume;
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=schedule)
    @test length(trace_before_resume) == length(trace_after_empty_resume)
    @test all(trace_before_resume .== trace_after_empty_resume)
    trace_after_empty_resume_2 = deepcopy(trace_before_resume)
    RepeaterPlacement.resume_simulated_annealing!(trace_after_empty_resume_2)
    @test length(trace_before_resume) == length(trace_after_empty_resume_2)
    @test all(trace_before_resume .== trace_after_empty_resume_2)
    trace_after_extension = deepcopy(trace_before_resume)
    RepeaterPlacement.resume_simulated_annealing!(trace_after_extension;
        num_epochs=num_epochs + 1)
    @test length(trace_after_extension) == length(trace_before_resume) + epoch_size
    for sol in trace_after_extension
        @test sol.X.kwargs[:num_epochs] == num_epochs + 1
        @test typeof(sol.X) == CostFunctionWithCacheAndKwargs
        @test 0. < sol.p[begin] < 11
        @test sol.X.cache[:cache_populated] == true
        @test sol.X.cache[:time] ≥ 0.
        @test sol.X.kwargs[:epoch_size] == epoch_size
        @test sol.X.kwargs[:schedule] == schedule
        @test sol.X.kwargs[:temperature] ≥ 0.
        @test 1 ≤ sol.X.kwargs[:epoch] ≤ num_epochs + 1
        @test sol.X.kwargs[:test_kwarg] == "test"
        @test sol.X.cache[:cache_from_evaluation] == true
    end
    # error for different epoch size
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume);
        num_epochs=num_epochs, epoch_size=epoch_size + 1, schedule=schedule)
    # error for different schedule
    @test_throws ArgumentError RepeaterPlacement.resume_simulated_annealing!(
        deepcopy(trace_before_resume);
        num_epochs=num_epochs, epoch_size=epoch_size, schedule=Exp(start=100., decay=0.8))
    # simulate an optimization stopping in the middle of an epoch
    trace_unfinished = deepcopy(trace_before_resume[1:end - floor(Int, epoch_size / 2)])
    trace_unfinished_resumed = deepcopy(trace_unfinished)
    RepeaterPlacement.resume_simulated_annealing!(trace_unfinished_resumed)
    @test length(trace_unfinished) < length(trace_unfinished_resumed)
    @test length(trace_unfinished_resumed) == length(trace_before_resume)
    @test trace_unfinished[end].p[begin] ≠ trace_unfinished_resumed[end].p[begin]

    # saving
    global saved
    saved = 0
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule, save=RepeaterPlacement.DoNotSave, savebase="test")
    @test saved == 0
    saved = 0
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule, save=RepeaterPlacement.SaveAtTheEnd, savebase="test")
    @test saved == 1
    saved = 0
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule, save=RepeaterPlacement.SaveEveryEpoch,
        savebase="test")
    @test saved == num_epochs
    saved = 0
    trace = RepeaterPlacement.simulated_annealing(cost_fct, initial_solution, num_epochs,
        epoch_size; schedule=schedule, save=RepeaterPlacement.SaveEveryIteration,
        savebase="test")
    @test saved == num_epochs * epoch_size
end
