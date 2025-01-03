
@testset "Creating prune fct" begin
    dominance_f(p1, p2) = path_length(p1) < path_length(p2)
    prune_f = create_prune_fct_from_dominance(dominance_f)
    path_1 = Path([1, 2, 3], [1., 1.])
    path_2 = Path([1, 2, 3], [2., 1.])
    path_3 = Path([1, 2, 3], [1., 2.])
    @test dominance_f(path_1, path_2) == true
    @test dominance_f(path_1, path_3) == true
    paths = [path_1, path_2, path_3]
    pruned_paths = prune_f(paths)
    @test pruned_paths == [path_1]
end

@testset "directly dominates" begin
    @test directly_dominates([1., 1., 1.], [2., 1., 2.])
    @test directly_dominates(Path([1, 2, 3, 4], [1., 1., 1.]),
        Path([5, 6, 7, 8], [2., 1., 2.]))
    @test ! directly_dominates([2., 2., 2.], [3., 1., 3.])
end

@testset "completely dominates" begin
    @test ! completely_dominates([1.0, 1.0, 1.0], [100.0])
    @test completely_dominates([1.0, 1.0, 1.0], [2.0, 1.0, 2.0])
    @test completely_dominates(Path([1, 2, 3, 4], [1.0, 1.0, 1.0]),
        Path([5, 6, 7, 8], [2.0, 1.0, 2.0]))
    @test ! completely_dominates([2.0, 2.0, 2.0], [3.0, 1.0, 3.0])
    @test completely_dominates([3.0, 2.0, 1.0], [4.0, 3.0, 2.0])
    @test completely_dominates([3.0, 2.0, 1.0], [2.0, 3.0, 4.0])
    @test completely_dominates([2.0, 2.0, 2.0], [1.0, 3.0, 3.0, 3.0, 1.0])
    @test ! completely_dominates([2.0, 2.0], [3.0, 1.0, 3.0])
end

@testset "incompletely dominates" begin
    @test ! incompletely_dominates([1.0, 1.0, 1.0], [100.0])
    @test incompletely_dominates([1.0, 1.0, 1.0], [2.0, 1.0, 2.0])
    @test completely_dominates(Path([1, 2, 3, 4], [1.0, 1.0, 1.0]),
        Path([5, 6, 7, 8], [2.0, 1.0, 2.0]))
    @test incompletely_dominates([1.0, 1.0, 1.0], [1.0, 2.0, 2.0, 2.0])
    @test incompletely_dominates([1.0, 1.0, 1.0], [2.0, 2.0, 2.0, 1.0])
    @test ! incompletely_dominates([2.0, 2.0, 2.0], [3.0, 1.0, 3.0])
    @test ! incompletely_dominates([2.0, 2.0, 2.0], [1.0, 3.0, 1.0, 3.0])
    @test ! incompletely_dominates([2.0, 2.0, 2.0], [3.0, 1.0, 3.0, 1.0])
    @test incompletely_dominates([3.0, 2.0, 1.0], [4.0, 3.0, 2.0])
    @test incompletely_dominates([3.0, 2.0, 1.0], [2.0, 3.0, 4.0])
    @test incompletely_dominates([3.0, 2.0, 1.0], [1.0, 4.0, 1.0, 3.0, 1.0, 2.0])
    @test incompletely_dominates([2.0, 2.0, 2.0], [1.0, 3.0, 3.0, 3.0, 1.0])
    @test incompletely_dominates([2.0, 2.0], [3.0, 1.0, 3.0])
    @test incompletely_dominates([2.0, 2.0], [1.0, 3.0, 1.0, 3.0, 1.0])
    @test ! completely_dominates([2.0, 2.0], [1.0, 3.0, 1.0, 3.0, 1.0])
end