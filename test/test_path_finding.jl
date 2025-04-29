@testset "Path" begin
    # make from nodes and lengths
    p = Path([1, 2, 3], [1.0, 2.0])
    @test p.nodes == nodes(p) == [1, 2, 3]
    @test p.lengths == [1.0, 2.0]
    @test edges(p) == [(1, 2), (2, 3)]
    @test length(p) == 2
    @test path_length(p) == 3.

    # make from graph
    g = SimpleWeightedGraph(4)
    add_edge!(g, 1, 2, 1.0)
    add_edge!(g, 2, 3, 2.0)
    add_edge!(g, 1, 3, 10.0)
    add_edge!(g, 3, 4, 15.0)
    p_from_graph = Path([1, 2, 3], g)
    @test p_from_graph.nodes == [1, 2, 3]
    @test p_from_graph.lengths == [1.0, 2.0]
    @test p_from_graph == p
    @test Path([1, 3], g).lengths == [10.0]
    @test_throws ArgumentError Path([1, 4], g)

    # make from coordinates
    coords = Coordinates()
    coords = add_end_node(coords, [0, 0])
    coords = add_end_node(coords, [10, 0])
    coords = add_repeater(coords, [5, 0])
    path = Path([1, 3, 2], coords)
    @test path == Path([1, 3, 2], [5.0, 5.0])

    # empty path
    p = Path()
    @test isempty(p)

    # adding
    p1 = Path([1, 2], [10.0])
    p2 = Path([2, 3], [20.0])
    @test p1 + p2 == Path([1, 2, 3], [10.0, 20.0])

    # equality
    p = Path([1, 2], [10.0])
    @test p == Path([1, 2], [10.0])
    @test p == Path([2, 1], [10.0])

    # hash
    p1 = Path([1, 2], [10.0])
    p2 = Path([2, 1], [10.0])
    @test hash(p1) == hash(p2)

    p1 = Path([1, 2, 1], [10.0, 20.0])
    p2 = Path([1, 2, 1], [20.0, 10.0])
    @test hash(p1) == hash(p2)
end

@testset "Path finding" begin
    # make a test graph that looks like a triangle:
    #      1
    #     / \
    #   1/   \3
    #   /     \
    #  2-------3
    #      1
    g_triangle = SimpleWeightedGraph(3)
    add_edge!(g_triangle, 1, 2, 1.)
    add_edge!(g_triangle, 2, 3, 1.)
    add_edge!(g_triangle, 1, 3, 3.)

    @testset "Path enumeration" begin
        # find all paths
        paths = RepeaterPlacement.find_all_paths(g_triangle, 1, 3)
        @assert paths == RepeaterPlacement.enumerate_paths(g_triangle, 1, 3)
        @assert length(paths) == 2
        p1, p2 = sort(paths, by=length)
        @test p1 == Path([1, 3], g_triangle)
        @test p2 == Path([1, 2, 3], g_triangle)

        # find k shortest paths
        paths = RepeaterPlacement.find_k_shortest_paths(g_triangle, 1, 3, 2)
        @assert paths == RepeaterPlacement.enumerate_paths(g_triangle, 1, 3, 2)
        @test length(paths) == 2
        p, = RepeaterPlacement.find_k_shortest_paths(g_triangle, 1, 3, 1)
        @test p == Path([1, 2, 3], g_triangle)
    end

    @testset "sample paths" begin
        paths = [Path([1, 2], [10.0]), Path([1, 2], [20.0])]
        path_cost_f(p; kwargs...) = path_length(p), Dict{Symbol,Any}(kwargs...)
        p, c, extra_cache = RepeaterPlacement.sample_path(paths, path_cost_f, x=5.)
        @test p == paths[1] && c == 10 && extra_cache[:x] == 5.
        p, c, extra_cache = RepeaterPlacement.sample_path(paths, path_cost_f,
            temperature=1.0)
        @test p isa Path
    end

    function test_minimize_path_length(f; kwargs...)
        h(p; k...) = path_length(p), Dict{Symbol, Any}(k...)
        p, cost, cache = f(h, g_triangle, 1, 3; y=8., kwargs...)
        @test p == Path([1, 2, 3], g_triangle)
        @test cost == path_length(p)
        @test cache isa Dict{Symbol, Any}
        @test cache[:y] == 8.
    end

    function test_minimize_one_over_lengths(f; kwargs...)
        h(p; k...) = sum(1 / l for l in p.lengths), Dict{Symbol, Any}(k...)
        p, cost, cache = f(h, g_triangle, 1, 3; z="z", kwargs...)
        @test p == Path([1, 3], g_triangle)
        @test cost == 1 / 3.
        @test cache isa Dict{Symbol, Any}
        @test cache[:z] == "z"
    end

    function test_non_isotonic(f; kwargs...)
        # non-isotonic but monotonic cost function
        # definition of isotonicity:
        # if cost(p1) < cost(p2), then cost(p1 + e) < cost(p2 + e)
        # here we engineer a counter example and see if we find the shortest path
        function cost_f(path::Path; kwargs...)
            diffs = [path.lengths[i + 1] - path.lengths[i] for i in 1:length(path) - 1]
            sum(diffs.^2), Dict{Symbol, Any}(:length => length(path))
        end
        # the cost is the square of the differences in lengths between adjacent edges
        # if the first edge has length a, the second x, and the third b,
        # then it is optimal to have x = (b - a) / 2, i.e., take even steps from a to b
        # however, the cost of the path formed by the first two edges is minimized by
        # x = a, i.e., by having the second edge as equal as possible to the first
        # this is an example demonstrating that the cost function is not isotonic
        
        g_iso = SimpleWeightedGraph(5)
        add_edge!(g_iso, 1, 2, 1.)
        add_edge!(g_iso, 2, 3, 2.)
        add_edge!(g_iso, 3, 4, 10.)
        add_edge!(g_iso, 1, 5, 1.)
        add_edge!(g_iso, 5, 3, 5)
        # this is the following graph:
        # 1 -- 2
        # |    |
        # 5 -- 3 -- 4

        # cost of 1-2-3 is smaller than cost of 1-5-3
        @test cost_f(Path([1, 2, 3], g_iso))[2] == Dict{Symbol, Any}(:length => 2)
        @test cost_f(Path([1, 2, 3], g_iso))[1] < cost_f(Path([1, 5, 3], g_iso))[1]
        # however, cost of 1-2-3-4 is larger than cost of 1-5-3-4
        @test cost_f(Path([1, 2, 3, 4], g_iso))[1] > cost_f(Path([1, 5, 3, 4], g_iso))[1]
        # test that it is still behaving monotonically though
        @test cost_f(Path([1, 2, 3], g_iso))[1] < cost_f(Path([1, 2, 3, 4], g_iso))[1]
        @test cost_f(Path([1, 5, 3], g_iso))[1] < cost_f(Path([1, 5, 3, 4], g_iso))[1]
        # test that we get the best path
        p, cost, cache = f(cost_f, g_iso, 1, 4; kwargs...)
        @test p == Path([1, 5, 3, 4], g_iso)
        @test cost == (10 - 5)^2 + (5 - 1)^2
        @test cache isa Dict{Symbol, Any}
        @test cache[:length] == length(p)
    end

    @testset "Path finding through enumeration" begin
        f = path_finding_through_enumeration
        test_minimize_path_length(f)
        prune_f = create_prune_fct_from_dominance(incompletely_dominates)
        test_minimize_path_length(f, prune_fct=prune_f)
        prune_f = create_prune_fct_from_dominance(completely_dominates)
        test_minimize_path_length(f, prune_fct=prune_f)
        # TODO should probabilistically fail for high temperature * boltzmann_constant
        # can we test whether it sometimes fails the test in that case?
        test_minimize_one_over_lengths(f)
    end

    @testset "Extended Dijkstra" begin
        # test retreiving path from "previous nodes" dict
        prev_nodes = Dict(
            3 => 2,
            2 => 1,
            1 => nothing
        )
        p = RepeaterPlacement._path_from_previous_nodes(prev_nodes, 3, g_triangle)
        @test p == Path([1, 2, 3], g_triangle)

        # test extended dijkstra when it reduces to regular dijkstra
        test_minimize_path_length(extended_dijkstra)

        # test extended dijkstra when minimizing sum(1/weight)
        test_minimize_one_over_lengths(extended_dijkstra)
    end

    @testset "monotonic-cost pathfinding" begin
        for method in [DepthFirstSearch(), BreadthFirstSearch(), WavefrontSearch()]
            f = (args...; kwargs...) -> path_finding_for_monotonic_costs(args...;
                method=method, kwargs...)
            test_minimize_path_length(f)
            # test if it also works if we already know the best path is below some value
            # there are two paths: one length 2, one length 3
            test_minimize_path_length(f, upper_bound=2.)
            test_minimize_path_length(f, upper_bound=2.5)
            test_minimize_path_length(f, upper_bound=3.)
            test_minimize_path_length(f, upper_bound=5.)
            test_minimize_path_length(f, upper_bound=3., prune_edges=true)
            # TODO can we test whether a test fails?
            # @test_throws test_minimize_path_length(f, upper_bound=1.)
            # we know that the best path is 1-2-3, so we can use that as guess
            test_minimize_path_length(f, node_sequence_guesses=[[1, 2, 3]])
            test_minimize_path_length(f, node_sequence_guesses=[[1, 2, 3]],
                prune_edges=true)
            # we can also give the worse path as guess, or both!
            test_minimize_path_length(f, node_sequence_guesses=[[1, 3]])
            test_minimize_path_length(f, node_sequence_guesses=[[1, 2, 3], [1, 3]])
            test_minimize_one_over_lengths(f)  # note: would fail if we prune edges
            test_non_isotonic(f)
            test_non_isotonic(f, prune_edges=true)
        end
    end
end

@testset "priority queue for different search methods" begin

    meth_to_type(::Type{DepthFirstSearch}) = Tuple{Int, Float64}
    meth_to_type(::Type{BreadthFirstSearch}) = Tuple{Int, Float64}
    meth_to_type(::Type{WavefrontSearch}) = Float64
    for M in [DepthFirstSearch, BreadthFirstSearch, WavefrontSearch]
        qs = [RepeaterPlacement.PathQueue{M}(), RepeaterPlacement.PathQueue(M())]
        for q in qs
            @test q isa RepeaterPlacement.PathQueue{M, meth_to_type(M)}
            @test q.priority_queue isa DataStructures.PriorityQueue{Path, meth_to_type(M)}
        end
    end

    p1 = Path([1, 2], [10.])
    c1 = 1.
    p2 = Path([1, 2, 3], [10., 10.])
    c2 = 2.
    p3 = Path([1, 3], [15.])
    c3 = 3.
    p4 = Path([1, 3, 4], [15., 15.])
    c4 = 4.

    function create_queue(method)
        q = RepeaterPlacement.PathQueue{method}()
        enqueue!(q, p1, c1)
        enqueue!(q, p2, c2)
        enqueue!(q, p3, c3)
        enqueue!(q, p4, c4)
        q
    end

    # Test DepthFirstSearch
    q = create_queue(DepthFirstSearch)
    @test dequeue_pair!(q) == (p2 => c2)  # deepest, lowest cost
    @test dequeue_pair!(q) == (p4 => c4)  # deepest, higher cost
    @test dequeue_pair!(q) == (p1 => c1)  # shallowest, lowest cost
    @test dequeue_pair!(q) == (p3 => c3)  # shallowest, higher cost
    q = create_queue(DepthFirstSearch)
    @test dequeue!(q) == p2
    @test dequeue!(q) == p4
    @test dequeue!(q) == p1
    @test dequeue!(q) == p3

    # Test BreadthFirstSearch
    q = create_queue(BreadthFirstSearch)
    @test dequeue_pair!(q) == (p1 => c1)  # shallowest, lowest cost
    @test dequeue_pair!(q) == (p3 => c3)  # shallowest, higher cost
    @test dequeue_pair!(q) == (p2 => c2)  # deepest, lowest cost
    @test dequeue_pair!(q) == (p4 => c4)  # deepest, higher cost
    q = create_queue(BreadthFirstSearch)
    @test dequeue!(q) == p1
    @test dequeue!(q) == p3
    @test dequeue!(q) == p2
    @test dequeue!(q) == p4

    # Test WavefrontSearch
    q = create_queue(WavefrontSearch)
    # just by order of cost
    @test dequeue_pair!(q) == (p1 => c1)
    @test dequeue_pair!(q) == (p2 => c2)
    @test dequeue_pair!(q) == (p3 => c3)
    @test dequeue_pair!(q) == (p4 => c4)
    q = create_queue(WavefrontSearch)
    @test dequeue!(q) == p1
    @test dequeue!(q) == p2
    @test dequeue!(q) == p3
    @test dequeue!(q) == p4
end

@testset "grow path" begin
    g = SimpleWeightedGraph(4)
    add_edge!(g, 1, 2, 1.0)
    add_edge!(g, 2, 3, 10.0)
    add_edge!(g, 2, 4, 100.0)
    path = Path([1, 2], g)
    paths = RepeaterPlacement.grow_path(path, g)
    expected_paths = [Path([1, 2, 3], g), Path([1, 2, 4], g)]
    @test (paths[1] == expected_paths[1] && paths[2] == expected_paths[2]) ||
        (paths[1] == expected_paths[2] && paths[2] == expected_paths[1])
end

@testset "prune graph" begin
    g = SimpleWeightedGraph(4)
    add_edge!(g, 1, 2, 1.0)
    add_edge!(g, 2, 3, 2.0)
    add_edge!(g, 3, 4, 3.0)
    add_edge!(g, 4, 5, 4.0)
    add_edge!(g, 5, 6, 5.0)
    cost_fct(p; k...) = path_length(p), Dict{Symbol, Any}()

    #extreme cases
    @test g == RepeaterPlacement.prune_graph(g, cost_fct, 10.)
    @test ne(RepeaterPlacement.prune_graph(g, cost_fct, 0.)) == 0

    # prune part of the edges
    g_pruned = RepeaterPlacement.prune_graph(g, cost_fct, 3.5)
    @test ne(g_pruned) == 3
    @test has_edge(g_pruned, 1, 2)
    @test has_edge(g_pruned, 2, 3)
    @test has_edge(g_pruned, 3, 4)
    @test g_pruned == RepeaterPlacement.prune_graph(g, cost_fct, 4.)
end

@testset "path from linear interpolation" begin
    f = RepeaterPlacement.path_from_linear_interpolation
    coords = Coordinates()
    coords = add_end_node(coords, [0, 0])
    coords = add_end_node(coords, [10, 0])
    coords = add_repeater(coords, [5, 0])
    path = f(coords, 1, 2, 0)
    @test path == Path([1, 2], [10.])
    path = f(coords, 1, 2, 1)
    @test path == Path([1, 3, 2], [5., 5.])
    path = f(coords, 1, 2)
    @test path == Path([1, 3, 2], [5., 5.])
    coords = Coordinates()
    coords = add_end_node(coords, [0, 0])
    coords = add_end_node(coords, [100, 0])
    for i in 1:99
        coords = add_repeater(coords, [i, 0])
    end
    path = f(coords, 1, 2, 9)  # 10 segments of 10 length
    @test path_length(path) == 100
    @test length(path) == 10
    for l in path.lengths
        @test l == 10
    end
end
