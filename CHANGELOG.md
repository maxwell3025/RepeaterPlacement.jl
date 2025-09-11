CHANGELOG
=========

Upcoming
--------
- Fixes in adding `Path`s: can be empty, errors actually thrown on mismatch.

2025-07-02 (v.0.4.0)
--------------------
- Added `waxmann_graph()` and `build_waxmann_graph()` for constructing Waxmann graphs with associated node coordinates.
- Joint plotting of graphs and coordinates: `plot_graph()` method with graph and `Coordinates`, `plot_node_locations()` can take graphs to determine which (special) edges to draw.

2025-06-09 (v.0.3.0)
--------------------
- Added an optional `rng` parameter for each of the randomized graph generation functions.
- Fixed an issue with docstrings in `Base.hash(p::Path)`.

2025-04-29 (v.0.2.0)
--------------------
- Implemented `Base.hash()` for the `Path` type.

2025-01-02 (v.0.1.0)
--------------------
- First public version, migrated from a private repository.
