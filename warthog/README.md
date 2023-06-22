## The Warthog Pathfinding Library 
**Author:** Daniel Harabor and contributors (see CONTRIB for a full list)  
**Contact:** daniel dot harabor at monash dot edu  
---

To compile: `make fast`  
To debug: `make dev`  
To profile: `make debug`  

By default we compile a small set of solver programs: `warthog`, `roadhog` and `mapf`. These can be found and executed from 
`./build/<mktarget>/bin` where `<mktarget>` is the name of the make target.

The full list of programs can be found in the `./programs/` directory. Currently 
they are as follows:

- `warthog.cpp`: for solving grid-based pathfinding problems
- `roadhog.cpp`: for solving road-network pathfinding problems
- `mapf.cpp`: for solving multi-agent pathfinding problems
- `fifo.cpp`: road-network pathfinding using a FIFO for i/o
- `make_cpd.cpp`: create a Compressed Path Database for a given input graph
- `grid2graph.cpp`: convert a gridmap to an xy-graph
- `dimacs2xy`: convert a DIMACS graph to an xy-graph
- `dimacs2metis`: convert a DIMACS graph to the input format of the METIS 
graph partitioning library.

Below we briefly describe the use of the `warthog` binary. For other programs 
refer to the inbuilt instructions that are printed on execution.  

**Looking for benchmarks to run?** These are kept in a separate [repository](https://bitbucket.org/shortestpathlab/benchmarks).

### Example: Solving grid-based pathfinding problems 

Once compiled, run bin/warthog for a list of command line parameters. 
A simple case is the following:

./bin/warthog --scen orz700d.map.scen --alg jps --checkopt

This invokes JPS on all instances in the scenario file orz700d.map.scen
The parameter --checkopt is optional and may be omitted. It simply
compares the length of a returned path with the optimal length specified
in the scenario file.

Various metrics are printed to the screen during execution. One line per instance.
Metrics are: nodes expanded, nodes generated (i.e. put on open), nodes touched
(i.e. evaluated, possibly resulting in an priority update), search time 
in microseconds (wallclock time) and path cost.

### Other options for the warthog program 

--alg [name]
Used to specify a named search algorithm.

--checkopt
Set this parameter to compare the length of each computed path against an
optimal length value specified by the scenario file at hand.

--gen [map file]
Used to generate random experiments over the specified map file.

--help
Set this parameter to print all available program options.

--scen [file]
Used to specify a scenario file for experiments.

--verbose
Set this parameter to print debugging information (use in conjunction with 
'make dev').

--wgm
Set this parameter to treat the map as a weighted-cost grid 
(cf. uniform cost). Under this model all tiles (save for explicit obstacles,
denoted by the input map using character '@') are considered traversable and 
have an associated cost equal to the ascii value used to describe the tile.
