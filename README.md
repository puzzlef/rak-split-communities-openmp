Design of OpenMP-based Parallel [Label Propagation Algorithm (LPA)] for [community detection], \
with no internally-disconnected communities.

Community detection is the problem of identifying tightly connected clusters of nodes within a network. Efficient parallel algorithms for this play a crucial role in various applications, especially as datasets expand to significant sizes. The Label Propagation Algorithm (LPA) is commonly employed for this purpose due to its ease of parallelization, rapid execution, and scalability. However, it may yield internally disconnected communities. This technical report introduces GSL-LPA, derived from our parallelization of LPA, namely GVE-LPA. Our experiments on a system with two 16-core Intel Xeon Gold 6226R processors show that GSL-LPA not only mitigates this issue but also surpasses FLPA, igraph LPA, and NetworKit LPA by `55x`, `10,300x`, and `5.8x`, respectively, achieving a processing rate of `844M` edges/s on a `3.8B` edge graph. Additionally, GSL-LPA scales at a rate of `1.6x` for every doubling of threads.

Below we plot the time taken by [FLPA], [igraph] LPA, [NetworKit] LPA, and GSL-LPA on 13 different graphs. GSL-LPA surpasses FLPA, igraph LPA, and NetworKit LPA by `55×`, `10,300×`, and `5.8×` respectively, achieving a processing rate of `844M` edges/s on a `3.8𝐵` edge graph.

[![](https://i.imgur.com/EmuwFSf.png)][sheets-o1]

Below we plot the speedup of GSL-LPA wrt FLPA, igraph LPA, and NetworKit LPA.

[![](https://i.imgur.com/dAgIdcx.png)][sheets-o1]

Next, we plot the modularity of communities identified by FLPA, igraph LPA, NetworKit LPA, and GSL-LPA. GSL-LPA on average obtains `7.1%` / `0.7%` higher modularity than FLPA and igraph LPA respectively (especially on road networks and protein k-mer graphs), and `3.6%` lower modularity than NetworKit LPA (especially on protein k-mer graphs).

[![](https://i.imgur.com/dqgA3ws.png)][sheets-o1]

Finally, we plot the fraction of disconnected communities identified by each implementation. Absence of bars indicates the absence of disconnected communities. GSL-LPA does not identify any communities that are internally-disconnected. However, on average, FLPA, igraph LPA, and NetworKit LPA exhibit fraction of disconnected communities amounting to `2.3×10^−3`, `1.2×10^−3`, and `1.8×10^−2`, particularly on web graphs and social networks.

[![](https://i.imgur.com/MG2HtXp.png)][sheets-o1]

Refer to our technical report for more details: [GSL-LPA: Fast Label Propagation Algorithm (LPA) for Community Detection with no Internally-Disconnected Communities][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.


[Label Propagation Algorithm (LPA)]: https://arxiv.org/abs/0709.2938
[FLPA]: https://github.com/vtraag/igraph/tree/flpa
[igraph]: https://github.com/igraph/igraph
[NetworKit]: https://github.com/networkit/networkit
[community detection]: https://en.wikipedia.org/wiki/Community_search
[sheets-o1]: https://docs.google.com/spreadsheets/d/18gPvelMwvANShEbE-zzC1Ui9wCdWja5xPS2jCVcBec0/edit?usp=sharing
[report]: https://arxiv.org/abs/2403.01261

<br>
<br>


### Code structure

The code structure of GSL-LPA is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/rak.hxx: LPA/RAK community detection algorithm functions
- inc/rakSplit.hxx: GSL-LPA community detection algorithm functions
- inc/split.hxx: Algorithms to split internally-disconnected communities
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetricize.hxx: Graph Symmetricization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


## References

- [Near linear time algorithm to detect community structures in large-scale networks; Raghavan et al. (2007)](https://arxiv.org/abs/0709.2938)
- [The University of Florida Sparse Matrix Collection; Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How to import VSCode keybindings into Visual Studio?](https://stackoverflow.com/a/62417446/1413259)
- [Configure X11 Forwarding with PuTTY and Xming](https://www.centlinux.com/2019/01/configure-x11-forwarding-putty-xming-windows.html)
- [Installing snap on CentOS](https://snapcraft.io/docs/installing-snap-on-centos)

<br>
<br>


[![](https://i.imgur.com/7QLfaW3.jpg)](https://www.youtube.com/watch?v=IwiYQILYXDQ)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/rak-split-communities-openmp)

[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
