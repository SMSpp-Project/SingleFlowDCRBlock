# Instances

The instances of the module, as netCDF files of the `Block` they describe,
which the tests of the module read. They are not in the repository: the
build downloads `nc4.tgz` from the Package Registry of the project, under the
version that `DATA_VERSION` in `CMakeLists.txt` names, and extracts it here
(target `extract_dcr_nc4`).

They come from 307 networks in four sets, each in a directory of its own,
since two sets may have networks of the same name (Abilene is both in sndlib
and in topo, in two different versions, and so are the ten GARR ones of garr
and of topo):

| set | networks | flows per network |
|---|---|---|
| `garr` | 10, the GARR network from 1999 to 2010 | 240 to 2970 |
| `sndlib` | 23, from SNDlib | 24 to 1482 |
| `topo` | 260, from the Internet Topology Zoo | 12 to 38612 |
| `waxman` | 14, random Waxman graphs | 664 to 39800 |

and for each network:

| directory | Block | what |
|---|---|---|
| `nc4/single/<set>` | `SingleFlowDCRBlock` | flow `<i>` = 0..9 on its own, `<network>_<i>.nc4`, each arc with its own capacity and cost |
| `nc4/multi/<set>` | `MultiFlowDCRBlock` | the first `<k>` = 1..5 flows, `<network>_<k>.nc4`, sharing the capacity of the arcs |

A formulation holding all the flows of a network is not what one solves
exactly, hence the multi-flow instances have the first few; the mutual
capacity of an arc depends on how many flows share it [see
`MultiFlowDCRBlock::load()`], so that an instance with its first k flows is a
problem of its own and not a part of the one with all of them, and with one
flow it is infeasible on most networks.

The textual instances they are written from are in `txt.tgz`, in the same
Package Registry under the same version (extracted into `txt/`, it is what
`make-nc4` reads): the multi-file format of `MultiFlowDCRBlock::load()`, with
the flows and the MTU in the `.dcr` (the "raw" one of `tools/dcr2nc4`), and,
for sndlib and waxman, the CSV of the experiments the sets come from, one
row per flow with the optimal value (`OVCplexInt`). `make-nc4` rebuilds
this directory out of them, and with `-a` it also writes the instances with
all the flows of each network, which are not distributed, each flow being a
group of the file.
The single-flow instances are the ones of the experiments the sets come
from, whose optimal values they reproduce.

New data go out as a new version, both archives together: raise
`DATA_VERSION`, then run `compress`, `upload-nc4`, `compress-txt` and
`upload-txt`.
