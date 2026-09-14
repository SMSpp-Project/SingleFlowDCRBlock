# SingleFlowDCRBlock

`SingleFlowDCRBlock` is a SMS++ module for Delay-Constrained Routing (DCR)
problems, i.e., the problem of routing one or more flows on a network so that
the worst-case end-to-end delay of each flow, as given by a network-calculus
formula, does not exceed its deadline, at minimum cost.

The module provides:

- `SingleFlowDCRBlock`, the :Block describing the DCR problem relative to a
  *single* flow: the binary routing variables, the reserved-rate and
  burst-delay variables, the flow conservation constraints, the delay
  constraint and the constraints linking rates to routing decisions. Two
  alternative formulations of the (nonlinear) burst-delay terms are
  supported, a Mixed-Integer Second-Order Cone one and a Perspective Cuts
  one in which the cones are dynamically outer-approximated by linear cuts;
  which of the two is constructed is selected by a `Configuration`;

- `MultiFlowDCRBlock`, the :Block describing the multicommodity version of
  the same problem: one `SingleFlowDCRBlock` sub-Block per flow, coupled by
  the mutual capacity constraints of the arcs;

- `DCRSolution`, the :Solution of `SingleFlowDCRBlock`;

- `SingleFlowDCRBendersSolver`, a :Solver for `SingleFlowDCRBlock`
  implementing a "Benders with nested Lagrangian relaxation" approach: the
  problem is decomposed along the minimum reserved rate, whose (scalar,
  nonconvex) value function is approximated by a cutting-plane scheme, each
  evaluation of which is a Lagrangian relaxation solved by a sequence of
  Shortest Path Tree computations.

The three "supporting" classes `DCR` (the abstract interface of a DCR
solver), `DCR_SPT` (two Shortest-Path-Tree heuristics for the single-flow,
single-path case) and `DCRLagrangianSolver` (the Lagrangian relaxation of
the delay constraint), together with the `SPT` Shortest Path Tree solver and
the `BenBound` cutting-plane scheme, are also part of the module: they do
not implement any SMS++ concept, they are the "physical" algorithms upon
which `SingleFlowDCRBendersSolver` is built.

See the documentation of each class for the details of the model and of the
algorithms.


## Getting started

These instructions will let you build the `SingleFlowDCRBlock` module on
your system.

The module also comes ready-made: `sudo apt install libsmspp-sfdcr-dev` from
the [PPA of the project](https://launchpad.net/~smspp-project/+archive/ubuntu/smspp),
which has `smspp-sfdcr` for its command-line tool as well, and `vcpkg install
"smspp[core,sfdcr]"` from the [SMS++ vcpkg
registry](https://gitlab.com/smspp/vcpkg-registry); `conda install -c
conda-forge smspp-project` and `brew install smspp`, from the [tap of the
project](https://github.com/SMSpp-Project/homebrew-smspp), carry the whole
framework. What follows is about building it yourself.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(SingleFlowDCRBlock)
target_link_libraries(<my_target> SMS++::SingleFlowDCRBlock)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `SingleFlowDCRBlock` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). These in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tests

The tests of the module live in the [tests](https://gitlab.com/smspp/tests)
project of the [Umbrella SMS++ Project](https://gitlab.com/smspp/smspp-project),
in the `SingleFlowDCRBlock` directory; see the `README.md` there.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new issue](https://gitlab.com/smspp/singleflowdcrblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Laura Galli**  
  Dipartimento di Informatica  
  Università di Pisa

- **Luca Mencarelli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Enrico Sorbera**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
