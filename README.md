Multi-threaded OpenMP-based Raghavan Albert Kumara ([RAK]) algorithm for
[community detection].

This is an implementation of a popular label-propagation based community
detection algorithm called **Raghavan Albert Kumara (RAK)**. Here, every node is
initialized with a unique label and at every step each node adopts the label
that most of its neighbors currently have. In this iterative process densely
connected groups of nodes form a consensus on a unique label to form
communities.

<br>
<br>


### Optimizations

#### OpenMP schedule

It appears `schedule(dynamic, 2048)` is the best choice.

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-schedule),
> [output](https://gist.github.com/wolfram77/5cc08eaceda9c523df9ac6e61ebd8f18), or
> [sheets][sheets-o1].

[![](https://i.imgur.com/2zJWRfx.png)][sheets-o1]
[![](https://i.imgur.com/nIaHUdY.png)][sheets-o1]

[sheets-o1]: https://docs.google.com/spreadsheets/d/1yIfh02VFVRU2_RMeh_mVS-djhAb5NA0N1dDeptCpC_w/edit?usp=sharing


#### Limiting number of iterations

It appears allowing a **maximum** of `20` **iterations** is ok.

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-max-iterations),
> [output](https://gist.github.com/wolfram77/00eec29e7c4f4debf6bc30a484a6d0a8), or
> [sheets][sheets-o2].

[![](https://i.imgur.com/eiH1VZz.png)][sheets-o2]
[![](https://i.imgur.com/m8sPfc6.png)][sheets-o2]

[sheets-o2]: https://docs.google.com/spreadsheets/d/1OpyFI0z89nsUFMPP8v9Dnr-ncZ0MMdDgxz7oFbr5PGc/edit?usp=sharing


#### Adjusting strictness

This controls how *equal-weighted labels* are treated. With **strict RAK**, we
always pick the *first label* with maximum total weight. With **non-strict**
**RAK**, we *randomly pick* among the top labels in case of a clash.

It appears **strict RAK** is faster (also better modularity).

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-strictness),
> [output](https://gist.github.com/wolfram77/aa67c7533739060e6ee9e588d83561ea), or
> [sheets][sheets-o3].

[![](https://i.imgur.com/bvOFkrP.png)][sheets-o3]
[![](https://i.imgur.com/PliysOi.png)][sheets-o3]

[sheets-o3]: https://docs.google.com/spreadsheets/d/1SFFb-RWlOSWV4_OE6-aHSG5FTHOEGSjoJfJ7l2YYzPA/edit?usp=sharing


#### Adjusting tolerance

It appears a **tolerance** of `0.05` is ok.

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-tolerance),
> [output](https://gist.github.com/wolfram77/a26f5ed0265a0a37709d86121ef7b1c0), or
> [sheets][sheets-o4].

[![](https://i.imgur.com/EdiNDJw.png)][sheets-o4]
[![](https://i.imgur.com/FnOvwga.png)][sheets-o4]

[sheets-o4]: https://docs.google.com/spreadsheets/d/1-F4zsNznUEVkjVIeH4rBdxwpUtUEHXbKm4B9qbsSjho/edit?usp=sharing


#### Vertex pruning

When a vertex changes its community, its marks its neighbors to be processed.
Once a vertex has been processed, it is marked as not to be processed. It helps
minimize unnecessary computation.

It appears with **vertex pruning** we get better timings.

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-pruning),
> [output](https://gist.github.com/wolfram77/a79bf240bf370ac2c96458418e774ed3), or
> [sheets][sheets-o5].

[![](https://i.imgur.com/WlQPFHl.png)][sheets-o5]
[![](https://i.imgur.com/PsePQGc.png)][sheets-o5]

[sheets-o5]: https://docs.google.com/spreadsheets/d/1Ze9qSB3mYTpboyqp5rj5NJCF-iY5dBEeR48EO_IiBlM/edit?usp=sharing


#### Hashtable design (data structure) for LPA iteration

One can use `P` `std::map`s (C++ inbuilt map) as the hastable for *LPA*. But
this has poor performance. So i use a **key-list** and a **full-size values**
**vector** (*collision-free*) we can dramatically improve performance. However if
the memory addresses of the hastables are **nearby**, even if each thread uses
its own hash table performance is not so high possibly due to false
cache-sharing (**Close-KV**). However if i ensure memory address are farther
away, the perf. improves (**Far-KV**).

It seems **Far-KV** has the best performance.

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/adjust-hashtable),
> [output](https://gist.github.com/wolfram77/8182cf5b270333aee8baa386a1198c61), or
> [sheets][sheets-o6].

[![](https://i.imgur.com/0YFz3CC.png)][sheets-o6]:
[![](https://i.imgur.com/rXuTxSf.png)][sheets-o6]

[sheets-o6]: https://docs.google.com/spreadsheets/d/1IA1N6Wh8DGuH4Y1zP73qhE0aBZyMlVLFfza9q1OGy00/edit?usp=sharing

<br>
<br>


### Main results

We combine the above *optimizations* and observe the performance of
**OpenMP-based RAK** on a number of graphs.

The input data used for this experiment is available from the
[SuiteSparse Matrix Collection]. This experiment was done with guidance
from [Prof. Kishore Kothapalli] and [Prof. Dip Sankar Banerjee].

> See
> [code](https://github.com/puzzlef/rak-communities-openmp/tree/input-large),
> [output](https://gist.github.com/wolfram77/5285c2a893f3ff8b1963eedee7be8967), or
> [sheets].

[![](https://i.imgur.com/3dvP2C5.png)][sheets]
[![](https://i.imgur.com/stkmVhD.png)][sheets]
[![](https://i.imgur.com/OrUoor8.png)][sheets]

[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[sheets]: https://docs.google.com/spreadsheets/d/1jkV3H6sBFgi6FTxJdvPnDtLPC-tp5kTeqtAfgS-ORT4/edit?usp=sharing

<br>
<br>


## References

- [Near linear time algorithm to detect community structures in large-scale networks; Usha Nandini Raghavan et al. (2007)](https://arxiv.org/abs/0709.2938)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How to import VSCode keybindings into Visual Studio?](https://stackoverflow.com/a/62417446/1413259)
- [Configure X11 Forwarding with PuTTY and Xming](https://www.centlinux.com/2019/01/configure-x11-forwarding-putty-xming-windows.html)
- [Installing snap on CentOS](https://snapcraft.io/docs/installing-snap-on-centos)

<br>
<br>


[![](https://i.imgur.com/7QLfaW3.jpg)](https://www.youtube.com/watch?v=IwiYQILYXDQ)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
