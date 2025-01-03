"""
    create_prune_fct_from_dominance(f)

From a function `f(p_1, p_2)` that indicates if `p_1` dominates `p_2`, create a function
that prunes a list of paths to only include the non-dominated paths.
"""
function create_prune_fct_from_dominance(f)
    function prune_fct(paths)
        for _ = 1:length(paths)
            p1 = popfirst!(paths)  # better to use `peel`?
            dominated = false
            for p2 in paths
                f(p2, p1) && (dominated = true; break)
            end
            dominated || push!(paths, p1)
        end
        paths
    end
end

"""
    directly_dominates(p1, p2)

`p2` is directly dominated by `p1` if all its segments are longer pairwise. Returns boolean.

Both paths should be represented as vectors of segments lengths or as `Path`s.

Direct domination is only defined for paths with the same number of segments;
this function will throw an error if `p1` and `p2` have a different number of segments.
By definition, `p1` directly dominates `p2` iff `p1[i] ≤ p2[i]` for all `i`.
If `p1` directly dominates `p2`, this suggests that the path `p1` will perform better
as a quantum-repeater chain than `p2`.
"""
function directly_dominates(p1, p2)
    length(p1) == length(p2) || throw(DimensionMismatch("p1 and p2 must have same length"))
    # We start with the assumption that p2 is dominated, and then try to disprove it.
    p2_dominated = true
    for (seg1, seg2) in zip(p1, p2)
        if seg1 > seg2
            p2_dominated = false
            break
        end
    end
    p2_dominated
end

directly_dominates(p1::Path, p2::Path) = directly_dominates(p1.lengths, p2.lengths)

"""
    completely_dominates(p1, p2)

Path `p1` completely dominates path `p2` if removing segments from the extremes of `p2` can
turn it into a path that is directly dominated by `p1`. Returns boolean.

Both paths should be represented as vectors of segments lengths or as `Path`s.

See also `directly_dominates` and `incompletely_dominates`.

Note that if `p1` completely dominates `p2`, it also incompletely dominates it.
Both relationships suggest that `p1` may be a better path for a quantum-repeater chain.
However, while this is not a strict implication for incomplete domination, it is for
complete domination.
"""
function completely_dominates(p1, p2)
    # If p1 has more segments than p2, p1 cannot dominate p2.
    # This is because p2 has not subpaths with the same number of segments as p1.
    if length(p1) > length(p2)
        return false
    elseif length(p1) == length(p2)
        subpaths = [p2, reverse(p2)]
    else
        diff = length(p2) - length(p1)
        subpaths = Vector{eltype(p2)}[]
        for i in 1:(diff + 1)
            subpath = p2[i:i+length(p1)-1]
            push!(subpaths, subpath, reverse(subpath))
        end
    end
    for subpath in subpaths
        directly_dominates(p1, subpath) && return true
    end
    false
end

completely_dominates(p1::Path, p2::Path) = completely_dominates(p1.lengths, p2.lengths)

"""
    incompletely_dominates(p1, p2)

Path `p1` incompletely dominates path `p2` if removing any set of segments from `p2` can
turn it into a path that is directly dominated by `p1`. Return boolean.

Both paths should be represented as vectors of segments lengths or as `Path`s.

See also `directly_dominates` and `completely_dominates`.

Note that if `p1` completely dominates `p2`, it also incompletely dominates it.
Both relationships suggest that `p1` may be a better path for a quantum-repeater chain.
However, while this is not a strict implication for incomplete domination, it is for
complete domination.

# Implementation
Sequentially, for each segment in `p1`, we look for the next segment in `p2` that is longer.
If we can do such a matching for each of the segments, that means that we could remove
all the segments in `p2` that were not matched with in order to create a path that is
directly dominated, and hence `p1` incompletely dominates `p2`.
If at some point we "run out" of segments in `p2`, this implies that it is impossible to
remove segments to create a path that is directly dominated by `p1`.
It remains possible though that segments can be removed to create a path that is directly
dominated by the reverse `p1`; this needs to be checked seperately.

"""
function incompletely_dominates(p1, p2)
    function _without_reversion(p1, p2)
        # This helper function tries to match segments between p1 and p2 but does not
        # account for the fact that p1 incompletely dominates p2 if this matching succeeds
        # for either p1 or the reverse of p1.
        i = 0
        for seg1 in p1
            while true
                i += 1
                i > length(p2) && return false  # we ran out of segments in p2
                seg1 ≤ p2[i] && break
            end
        end
        true
    end
    if length(p1) > length(p2)
        return false
    elseif length(p1) == length(p2)
        directly_dominates(p1, p2) || directly_dominates(reverse(p1), p2)
    else
        _without_reversion(p1, p2) || _without_reversion(reverse(p1), p2)
    end
end

incompletely_dominates(p1::Path, p2::Path) = incompletely_dominates(p1.lengths, p2.lengths)