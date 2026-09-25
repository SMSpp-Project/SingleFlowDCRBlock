# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the instances of the module, as netCDF files downloaded from the Package
  Registry into `data/nc4` [see `data/README.md`]: 140 `SingleFlowDCRBlock`,
  one per flow of 14 networks (the ten GARR ones, Abilene, Cogentco, Colt and
  w1_100_04), and 70 `MultiFlowDCRBlock`, each network with its first 1 to 5
  flows, which is what a formulation holding all of them solves exactly

- `tools/dcr2nc4`, which writes a `MultiFlowDCRBlock` out of the multi-file
  textual format that `load()` reads, or out of the raw one of the data sets
  it comes from, possibly keeping its first flows only, working in a
  temporary directory, since `load()` writes its own files in the current one

### Changed

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- `chg_costs()`, `chg_ucaps()` and `chg_dfcts()` take their data as a
  `std::span< const double >`, whose length they check against the Range or
  the Subset instead of reading past the end, and the first two are
  registered in the methods factory in that form too; the forms taking an
  iterator stay, and defer to the span ones

### Fixed

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- `is_sol_feasible()` declared infeasible a routing where an arc the flow
  does not use reserves a rate of the order of the tolerance of the Solver:
  r <= U x was checked against a tolerance relative to U x, i.e., an
  absolute one when x = 0, and it is now relative to U, the size of the
  terms of the constraint

- a `MultiFlowDCRBlock` could not be read back from the netCDF file it
  wrote: `serialize()` wrote the topology and the costs alone, leaving out
  the data of the flows, and `deserialize()` built no `SingleFlowDCRBlock`
  and kept the previous instance, so that generating the abstract
  representation of what it read crashed. The file now holds the mutual
  capacities, the data of each flow and the `SingleFlowDCRBlock` of each flow
  in a group of its own, which `deserialize()` rebuilds; the matrices of the
  costs, of the capacities and of the deficits of an earlier file are still
  read, and no longer written, the costs being in the flows

- `map_forward_Modification()` passed a change of the deficits on to the
  other `SingleFlowDCRBlock` as a change of its capacities, and, when this
  one had no capacities or no deficits, read them out of the empty vector
  rather than passing on the infinite capacities and the zero deficits it
  has

- `BenBound.cpp` no longer includes `unistd.h`, which it does not use and
  which MSVC does not have

## [0.1.1] - 2026-09-13

### Fixed

- `is_sol_feasible()` reads the DCRSolution, and the Variable are left alone

## [0.1.0] - 2026-09-12

### Added

- the standard SMS++ module build system: the CMake project (with the
  package configuration files and the version derived from the git tag),
  the CI for both GitLab and GitHub, and the `makefile` / `makefile-c` /
  `makefile-s` triple of every other module

- `SingleFlowDCRBlock::is_sol_feasible()`, which checks the solution held
  by a `DCRSolution` without touching the `Block`, i.e., without requiring
  its abstract representation to exist at all

- the accessors and the setters of `DCRSolution`, so that a `Solver` can
  fill it directly out of its own data structures

- `SingleFlowDCRBendersSolver::get_Solution()`, which does exactly that:
  the `DCRSolution` is filled out of what BenBound found, and the `Variable`
  of the `Block` are not written at all

- `SingleFlowDCRBlock::link_feasible()` and `delay_feasible()`, the two
  halves of the feasibility check that were nowhere to be found: the
  constraints tying the reserved rates to the routing decisions, and the
  end-to-end delay constraint

### Changed

- `SingleFlowDCRBlock::is_feasible()` checks all the constraints of the
  problem, the bounds of the variables and the two cone constraints
  comprised, and its "physical" branch reads the solution out of the
  `Variable` and checks it against the data, where it used to do nothing

- each of the four feasibility checks has a version taking the solution to
  be checked from the outside, which is what the "physical" version and
  `is_sol_feasible()` both boil down to

- `SingleFlowDCRBendersSolver::get_var_solution()` also writes the two
  "aggregate" variables `r_min` and `theta_min`, which it used to compute
  and drop; `is_DCR_feasible()` asks `SingleFlowDCRBlock::delay_feasible()`
  rather than repeating the network calculus formula, and the data of the
  `Block` are translated for BenBound in one place instead of two

- `SingleFlowDCRBendersSolver::has_var_solution()` says no unless the point
  BenBound stopped at is feasible for the DCR problem, which is now checked
  in full (the routing, the rates and the delay) and with a tolerance of
  1e-6 rather than the 1e-2 of the delay alone; `get_ub()` accordingly
  returns +INF where there is no feasible solution, in place of the
  (possibly optimistic) relaxed value of `BenBound::getHeurVal()`: neither
  that value nor the value of an infeasible point bounds the optimum

- `SingleFlowDCRBendersSolver::get_var_value()` computes the value of the
  solution BenBound found rather than of whatever the `Variable` of the
  `Block` happen to hold, which is what anybody could have written there

- `MultiFlowDCRBlock::is_feasible()` asks each commodity sub-`Block` about
  the commodity it describes, where it only checked the mutual capacity
  constraints of its own

- which of the two formulations is in use is asked for in the
  `Configuration` of the *static Constraint* of the `BlockConfig`, the one
  of the static `Variable` being still honoured as it used to be; the two
  methods that ask no longer disagree on the default, which is the SOCP
  formulation for both

- `SingleFlowDCRBlock::load()` takes the starting nodes before the ending
  ones, as its own documentation says

- `get_NodeDelays()` and `get_LinkDelays()` return a reference, as every
  other datum of the `Block` does

- the timer of the DCR solvers is the `std::chrono`-based `DCRtimer`,
  which builds everywhere, in place of the POSIX-only `OPTtimers`

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Removed

- the local copies of `OPTUtils.h`, `OPTtypes.h` and `OPTvect.h`, some
  3700 lines of somebody else's code that the module carried around

- `SingleFlowDCRBlock::is_feasible_flow()`, which did exactly the same
  checks as `is_feasible()` under a name saying otherwise

### Fixed

- the aggregate cone constraint was added to the abstract representation
  even when the "P/C" formulation was asked for, so that formulation was
  not linear at all and a `:MILPSolver` that does not do conic constraints
  could not solve it

- the debug output that the solvers wrote on `std::cout` at every
  iteration, and the `exit( 1 )` calls with which they killed the whole
  process where they should throw

- `is_feasible_instance()` read the node deficits, the arc capacities and
  the delays out of vectors that are allowed to be empty, and left the
  source and the sink uninitialised when the instance had none

- `DCR_SPT::DCRgetLink()` fell off the end of the function when the arc was
  not there, and `DCRLagrangianSolver` compared a `double` through the
  integer `abs()`

- `BenBound` moved a trial value of `r_min` that it had already solved for
  by shrinking it multiplicatively, which walks past the left endpoint of
  the window where an `r_min` can be feasible at all: it could therefore
  reserve less than the sustained rate of the flow, and hand back a point
  that is not a solution together with a value below the optimum. The
  point is now moved towards that endpoint, and stops there

[Unreleased]: https://gitlab.com/smspp/singleflowdcrblock/-/compare/0.1.1...develop
[0.1.1]: https://gitlab.com/smspp/singleflowdcrblock/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/singleflowdcrblock/-/tags/0.1.0
