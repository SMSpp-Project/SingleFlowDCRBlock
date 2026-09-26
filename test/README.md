# test

A tester for the `SingleFlowDCRBlock` module that needs nothing but the core
SMS++ library.

All its instances are tiny graphs built in memory, whose simple s-t paths are
few enough to be enumerated, so that every answer of the solvers of the module
is compared with the one the enumeration gives. `SPT` is compared with the
shortest of the paths, in the path it returns and in the label of every node,
on random graphs with non-negative costs and on random acyclic ones with
negative costs. `DCR_SPT` is compared with the paths that meet the deadline:
ERA-I has to say that an instance is infeasible exactly when no path meets it
with the whole capacity reserved on every arc, and what either heuristic
returns has to be a path meeting the deadline at the rates it reports, whose
cost is not below the continuous optimum (computed path by path, by bisection
on the multiplier of the delay constraint). `DCRLagrangianSolver`, for a given
r_min, has to give a Lagrangian value not above the optimum of the problem it
relaxes, and a heuristic value, with the rates that attain it, not below the
optimum of the DCR problem on the same reduced graph.

`SingleFlowDCRBlock` is checked in the feasibility of given solutions, the
default data of `load()`, the netCDF round trip (the file is written in the
working directory and removed), `is_feasible_instance()` against the
enumeration, and the physical `Modification` issued by the changes of costs,
capacities, source and sink and by closing and opening arcs, which a `Solver`
attached to the `Block` sees when it solves again; the delays and the deadline
cannot be changed but by loading the instance again.

The edge cases are checked on their own: source equal to sink, no path at all,
every path violating the deadline, a deadline met with equality, ties between
paths, zero-cost arcs, parallel arcs and closed arcs. The random instances are
drawn with fixed seeds.

The Benders solver of the module is not tested here, since it takes the
`Solver` of its master problem from a configuration and that needs a MILP
solver.

The exit code is the number of failed checks, and `All tests passed!!` is
printed when there is none. The `makefile` builds the executable including the
`SingleFlowDCRBlock` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
