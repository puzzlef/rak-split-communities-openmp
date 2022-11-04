Multi-threaded OpenMP-based Raghavan Albert Kumara ([RAK]) algorithm for
[community detection].

This is an implementation of a popular label-propagation based community
detection algorithm called **Raghavan Albert Kumara (RAK)**. Here, every node is
initialized with a unique label and at every step each node adopts the label
that most of its neighbors currently have. In this iterative process densely
connected groups of nodes form a consensus on a unique label to form
communities.

When there exist multiple communities with max weight, we randomly pick one of
them (**non-strict** approach), or pick only the first of them (**strict** approach).
The algorithm converges when `n%` of vertices dont change their community
membership (*tolerance*).

We continue with OpenMP implementation of RAK algorithm for community detection.
Each thread is given a *separate* hashtable, which it can use for choosing the
most popular label among its neighbors (by weight). I initially packed the
hashtables (for each thread) contiguously on a vector. However, i observed that
*allocating them separately* given almost *2x* performance (by using pointer to
vectors). OpenMP schedule is `auto` now, we can later choose the best if we
need.

[![](https://i.imgur.com/zLLrbnj.png)][sheetp]

<br>
<br>

Similar to [previous experiment], we adjust *tolerance* from `0.1` to `0.0001` and
compare the sequential and OpenMP implementations (non-strict, strict
approaches).

[![](https://i.imgur.com/10vwdJf.png)][sheetp]

[![](https://i.imgur.com/uhgVh1A.png)][sheetp]

[![](https://i.imgur.com/hH3sKyS.png)][sheetp]

<br>
<br>

On average, **OpenMP-based strict RAK appears to be better**, both in terms of
time and modularity (2X with 1->). Also we generally see it do better than
*sequential* approaches (possibly due to better randomization). For a
*tolerance* of `0.05`, *OpenMP-based strict RAK* (using *12* threads) is *6.75x*
faster than *sequential non-strict RAK*.

However, OpenMP implementations do not achieve better quality for `coAuthors*`
graphs. For *social networks*, OpenMP-based non-strict RAK achieves better
modularity than the strict version. Please see [previous experiment] for difference
between strict and non-strict versions. Given below is modularity on
`soc-LiveJournal1` graph.

[![](https://i.imgur.com/4GMolxZ.png)][sheetp]

<br>
<br>

As brefore, it seems to me a *tolerance* of `0.05` would be a good choice.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].


[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
[previous experiment]: https://github.com/puzzlef/rak-communities-seq
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 3985272 [directed] {} (symmetricize)
# OMP_NUM_THREADS=12
# [-0.000497 modularity] noop
# [00150.956 ms; 0003 iters.; 0.872759998 modularity] rakSeqStatic       {tolerance=1e-01}
# [00152.820 ms; 0003 iters.; 0.847129285 modularity] rakSeqStaticStrict {tolerance=1e-01}
# [00021.234 ms; 0003 iters.; 0.872025967 modularity] rakOmpStatic       {tolerance=1e-01}
# [00021.151 ms; 0003 iters.; 0.847373843 modularity] rakOmpStaticStrict {tolerance=1e-01}
# [00149.791 ms; 0003 iters.; 0.872759998 modularity] rakSeqStatic       {tolerance=5e-02}
# [00147.402 ms; 0003 iters.; 0.847129285 modularity] rakSeqStaticStrict {tolerance=5e-02}
# [00020.572 ms; 0003 iters.; 0.871711493 modularity] rakOmpStatic       {tolerance=5e-02}
# [00020.510 ms; 0003 iters.; 0.847401798 modularity] rakOmpStaticStrict {tolerance=5e-02}
# [00236.676 ms; 0005 iters.; 0.876379907 modularity] rakSeqStatic       {tolerance=1e-02}
# [00234.186 ms; 0005 iters.; 0.853610516 modularity] rakSeqStaticStrict {tolerance=1e-02}
# [00030.638 ms; 0005 iters.; 0.875127017 modularity] rakOmpStatic       {tolerance=1e-02}
# [00030.697 ms; 0005 iters.; 0.852154553 modularity] rakOmpStaticStrict {tolerance=1e-02}
# [00284.453 ms; 0006 iters.; 0.877100408 modularity] rakSeqStatic       {tolerance=5e-03}
# [00277.874 ms; 0006 iters.; 0.854387999 modularity] rakSeqStaticStrict {tolerance=5e-03}
# [00035.626 ms; 0006 iters.; 0.876359046 modularity] rakOmpStatic       {tolerance=5e-03}
# [00035.454 ms; 0006 iters.; 0.853212833 modularity] rakOmpStaticStrict {tolerance=5e-03}
# [00366.719 ms; 0008 iters.; 0.877606571 modularity] rakSeqStatic       {tolerance=1e-03}
# [00407.965 ms; 0009 iters.; 0.855403066 modularity] rakSeqStaticStrict {tolerance=1e-03}
# [00049.822 ms; 0009 iters.; 0.876784205 modularity] rakOmpStatic       {tolerance=1e-03}
# [00055.832 ms; 0010 iters.; 0.853845000 modularity] rakOmpStaticStrict {tolerance=1e-03}
# [00448.478 ms; 0010 iters.; 0.877725899 modularity] rakSeqStatic       {tolerance=5e-04}
# [00575.903 ms; 0013 iters.; 0.856024683 modularity] rakSeqStaticStrict {tolerance=5e-04}
# [00064.744 ms; 0011 iters.; 0.876972854 modularity] rakOmpStatic       {tolerance=5e-04}
# [00069.341 ms; 0013 iters.; 0.854331434 modularity] rakOmpStaticStrict {tolerance=5e-04}
# [00616.509 ms; 0014 iters.; 0.877830148 modularity] rakSeqStatic       {tolerance=1e-04}
# [00742.345 ms; 0017 iters.; 0.856140256 modularity] rakSeqStaticStrict {tolerance=1e-04}
# [00102.798 ms; 0020 iters.; 0.877253771 modularity] rakOmpStatic       {tolerance=1e-04}
# [00098.093 ms; 0019 iters.; 0.854501545 modularity] rakOmpStaticStrict {tolerance=1e-04}
#
# ...
```

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


[gist]: https://gist.github.com/wolfram77/d4503226d989c2752210df65ea12ec4d
[charts]: https://imgur.com/a/cYzo2Ai
[sheets]: https://docs.google.com/spreadsheets/d/1D7EpBMmnGlJlk0uUEqUTYWm5Gxc-AXVhSfxbDA56y8Y/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRLy5tronSINq10-myRK8M7ykPKwXF0AwwvssViiqbu3va6USoVncYppn6RzvxNqGw8ev2gIDQ1G7wA/pubhtml
