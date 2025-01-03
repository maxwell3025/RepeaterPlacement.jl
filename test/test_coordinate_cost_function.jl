end_node_coordinate_mat = [[0.0, 0.0] [10.0, 0.0]]
function f(x::Coordinates; index, kwargs...)
    endnode = end_node(x, index)
    rep = repeater(x, 1)
    cost = RepeaterPlacement.distance(endnode, rep)
    cache = Dict(:end_node_x_coord => endnode[1])
    :cost_function_value in keys(kwargs) &&
        (cache[:previous_cost] = kwargs[:cost_function_value])
    cost, cache
end
cost_fct_1 = CoordinateCostFunction(f, end_node_coordinate_mat; index=1)
end_node_coords = Coordinates()
end_node_coords = add_end_node(end_node_coords, [0.0, 0.0])
end_node_coords = add_end_node(end_node_coords, [10.0, 0.0])
end_node_coords = add_repeater(end_node_coords, [100.0, 100.0])
cost_fct_2 = CoordinateCostFunction(f, end_node_coords; index=1)
cost_fct_3 = copy(cost_fct_1)
rep_coord_mat = [[80.0, 0.0];;]
for cost_fct in [cost_fct_1, cost_fct_2, cost_fct_3]
    @test cost_fct.end_node_coordinate_matrix == end_node_coordinate_mat
    @test cost_fct.kwargs == Dict(:index => 1)
    cost = cost_fct(rep_coord_mat)
    @test cost == 80.0
    @test cost_fct.cache[:cost_function_value] == cost
    @test cost_fct.cache[:end_node_x_coord] == 0.0
    cost_fct.kwargs[:index] = 2
    cost = cost_fct(rep_coord_mat)
    @test cost == 70.0
    @test cost_fct.cache[:cost_function_value] == cost
    @test cost_fct.cache[:end_node_x_coord] == 10.0
    @test cost_fct.cache[:previous_cost] == 80.0
end
# also test if the cache is copied correctly (above doesn't test because eval populates it)
@test copy(cost_fct_1).cache == cost_fct_1.cache