Single-threaded CPU-based Raghavan Albert Kumara ([RAK]) algorithm for
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

[![](https://i.imgur.com/g6Zn58k.png)][sheetp]

<br>
<br>

In this experiment we adjust *tolerance* from `0.1` to `0.0001`.

[![](https://i.imgur.com/gv2JJPn.png)][sheetp]

[![](https://i.imgur.com/smaaCf0.png)][sheetp]

[![](https://i.imgur.com/XzQmY14.png)][sheetp]

<br>
<br>

Non-strict approach generally achieves better modularity, with certain
exceptions (`coAuthors*`, road networks). I think this has to do with graph
structure. See for example on `coAuthorsDBLP` graph.

[![](https://i.imgur.com/YrmCLil.png)][sheetp]

<br>
<br>

It seems to me a *tolerance* of `0.05` would be a good choice.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].


[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
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
# [-0.000497 modularity] noop
# [00157.008 ms; 0003 iters.; 0.872759998 modularity] rakSeqStatic       {tolerance=1e-01}
# [00157.307 ms; 0003 iters.; 0.847129285 modularity] rakSeqStaticStrict {tolerance=1e-01}
# [00156.578 ms; 0003 iters.; 0.872759998 modularity] rakSeqStatic       {tolerance=5e-02}
# [00152.985 ms; 0003 iters.; 0.847129285 modularity] rakSeqStaticStrict {tolerance=5e-02}
# [00246.655 ms; 0005 iters.; 0.876379907 modularity] rakSeqStatic       {tolerance=1e-02}
# [00246.319 ms; 0005 iters.; 0.853610516 modularity] rakSeqStaticStrict {tolerance=1e-02}
# [00290.881 ms; 0006 iters.; 0.877100408 modularity] rakSeqStatic       {tolerance=5e-03}
# [00286.677 ms; 0006 iters.; 0.854387999 modularity] rakSeqStaticStrict {tolerance=5e-03}
# [00384.414 ms; 0008 iters.; 0.877606571 modularity] rakSeqStatic       {tolerance=1e-03}
# [00420.077 ms; 0009 iters.; 0.855403066 modularity] rakSeqStaticStrict {tolerance=1e-03}
# [00469.554 ms; 0010 iters.; 0.877725899 modularity] rakSeqStatic       {tolerance=5e-04}
# [00602.378 ms; 0013 iters.; 0.856024683 modularity] rakSeqStaticStrict {tolerance=5e-04}
# [00648.129 ms; 0014 iters.; 0.877830148 modularity] rakSeqStatic       {tolerance=1e-04}
# [00776.319 ms; 0017 iters.; 0.856140256 modularity] rakSeqStaticStrict {tolerance=1e-04}
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


[![](https://i.imgur.com/oIHCg7z.jpg)](https://www.youtube.com/watch?v=N68Lha2Wa6U)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)


[gist]: https://gist.github.com/wolfram77/0ecfd0796a5fc8ee1c42bcb77b696aec
[charts]: https://imgur.com/a/YHySiIn
[sheets]: https://docs.google.com/spreadsheets/d/1UR7ZCDYA6Ed7Yi66IVpAQujChAiKOUbz2YH6urzQWhw/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vS8-QCUgT6zNRygu6hYNt2rPU2cZXvFr3Mr31OPa3bLOEJ2mEzbwSBOrI-DyyOML_Lc6FiPkbneL-yk/pubhtml
