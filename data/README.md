# Instances

The instances of the module, as netCDF files of the `Block` they describe,
which the tests of the module read. They are not in the repository: the
build downloads `nc4.tgz` from the Package Registry of the project, under the
version that `DATA_VERSION` in `CMakeLists.txt` names, and extracts it here
(target `extract_dcr_nc4`).

| directory | Block | what |
|---|---|---|
| `nc4/single` | `SingleFlowDCRBlock` | flow `<i>` = 0..9 of each network, `<network>_<i>.nc4`, 140 instances |
| `nc4/multi` | `MultiFlowDCRBlock` | the first `<k>` = 1..5 flows of each network, `<network>_<k>.nc4`, 70 instances |

The networks are 14: the ten GARR ones, from `Garr199901` to `Garr201001`,
and `Abilene`, `Cogentco`, `Colt` and `w1_100_04`. A network has from 110
(`Abilene`) to 38612 (`Cogentco`) flows, hence a formulation holding all of
them is not what one solves exactly, and the instances with the first few
flows are; the mutual capacity of an arc depends on how many flows share it
[see `MultiFlowDCRBlock::load()`], so that an instance with its first k
flows is a problem of its own and not a part of the one with all of them.

The multi-flow instances are written by `tools/dcr2nc4` out of the textual
format that `MultiFlowDCRBlock::load()` reads; `make-nc4` rebuilds the whole
directory out of the textual instances, which are kept in the tag
`archive/dcr-flow` of the `tests` repository, and with `-a` it also writes
the instances with all the flows of each network, which are not distributed:
each flow is a group of the file, and they are up to 38612. New data go out as a new
version: raise `DATA_VERSION`, then run `compress` and `upload-nc4`.
