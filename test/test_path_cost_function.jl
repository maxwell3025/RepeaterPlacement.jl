@testset "largest pathwise cost" begin
    @testset "multiple end-node pairs" begin
        # example with multiple end-node pairs
        coords = Coordinates()
        coords = add_end_node(coords, [0., 0.])
        coords = add_end_node(coords, [5., 0.])
        coords = add_end_node(coords, [10., 0.])
        # coords = add_repeater(coords, [5., 0.])
        # all the possible paths, one per pair of end nodes
        p12 = Path([1, 2], [5.])
        p13 = Path([1, 3], [10.])  # this is the worst one
        p23 = Path([2, 3], [5.])
        # p2 = Path([1, 3, 2], [5., 5.])
        # f_path(p::Path) = sum(1 / l for l in p.lengths)
        f_path(p::Path; kwargs...) = path_length(p), Dict{Symbol, Any}()
        cost, cache = RepeaterPlacement.largest_pathwise_cost(coords, f_path)
        @test cost == 10.
        @test cache[:path] == p13

        # use this example to also test the wrapper LargestPathwiseCostFct
        f = RepeaterPlacement.LargestPathwiseCostFct(f_path)
        @test path_cost_fct(f) == f_path
        cost2, cache2 = f(coords)
        @test cost2 == cost; @test cache2 == cache
        # and we also wrap it in a CoordinateCostFunction
        coord_cost_fct = CoordinateCostFunction(f, coords)
        @test path_cost_fct(coord_cost_fct) == f_path
        cost3 = coord_cost_fct(coords.repeater_coordinate_matrix)
        @test cost3 == cost
        cache[:cost_function_value] = cost
        @test coord_cost_fct.cache == cache
    end

    @testset "multiple paths" begin
        # example with one pair of end nodes but multiple paths
        coords = Coordinates()
        coords = add_end_node(coords, [0., 0.])
        coords = add_end_node(coords, [10., 0.])
        coords = add_repeater(coords, [5., 0.])
        f_path(p::Path; kwargs...) = sum(l^2 for l in p.lengths),
            Dict{Symbol, Any}(:num_reps => length(p) - 1)
        p_direct = Path([1, 2], [10.])  # cost = 100
        p_repeater = Path([1, 3, 2], [5., 5.])  # cost = 50, this is the best path
        cost, cache = RepeaterPlacement.largest_pathwise_cost(coords, f_path)
        @test cost == 50.; @test cache[:num_reps] == 1
    end

    @testset "use one specific path only" begin
        coords = Coordinates()
        coords = add_end_node(coords, [0., 0.]) # 1
        coords = add_end_node(coords, [10., 0.]) # 2
        coords = add_end_node(coords, [20., 0.])  # 3
        coords = add_end_node(coords, [10., 50.])  # 4
        coords = add_repeater(coords, [0., 5.])  # 5
        coords = add_repeater(coords, [10., 5.])  # 6
        coords = add_repeater(coords, [20., 5.])  # 7
        # many paths are possible, but we should only be using one
        f_path(p::Path; kwargs...) = path_length(p), Dict{Symbol, Any}()
        path = Path([1, 5, 7, 3], [1., 1., 1.])
        # this paths total length is 5 + 20 + 5 = 30
        # the segment lengths are wrong in the path specification,
        # should be updated automatically in the cost function call
        cost, cache = RepeaterPlacement.largest_pathwise_cost(coords, f_path,
            path=path, reuse_path=true)
        @test cost == 30.
        @test cache[:path] == Path([1, 5, 7, 3], [5., 20., 5.])
    end
end

@testset "mean pathwise cost" begin
    f_path(p::Path; kwargs...) = path_length(p), Dict{Symbol, Any}()

    # three pairs of end nodes
    coords = Coordinates()
    coords = add_end_node(coords, [0., 0.])
    coords = add_end_node(coords, [5., 0.])
    coords = add_end_node(coords, [10., 0.])
    expected_cost = (10. + 5. + 5.) / 3
    cost, cache = RepeaterPlacement.mean_pathwise_cost(coords, f_path)
    @test cost == expected_cost

    # test also wrappers
    f = RepeaterPlacement.MeanPathwiseCostFct(f_path)
    @test path_cost_fct(f) == f_path
    cost2, cache2 = f(coords)
    @test cost2 == cost; @test cache2 == cache
    # and we also wrap it in a CoordinateCostFunction
    coord_cost_fct = CoordinateCostFunction(f, coords)
    @test path_cost_fct(coord_cost_fct) == f_path
    cost3 = coord_cost_fct(coords.repeater_coordinate_matrix)
    @test cost3 == cost
    cache[:cost_function_value] = cost
    @test coord_cost_fct.cache == cache

    # one pair, two paths
    coords = Coordinates()
    coords = add_end_node(coords, [0., 0.])
    coords = add_end_node(coords, [1., 1.])
    coords = add_repeater(coords, [1., 0.])
    # direct path has length sqrt(2) (Pythagoras), via repeater length 2
    # we take the mean over all the best paths of each pair,
    # now there is only one pair so we just pick the best path
    cost, cache = RepeaterPlacement.mean_pathwise_cost(coords, f_path)
    @test cost == sqrt(2)
end