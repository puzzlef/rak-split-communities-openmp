#pragma once
#include <cstdint>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif
#include "dfs.hxx"
#include "bfs.hxx"

using std::vector;




#pragma region SPLIT DISCONNECTED COMMUNITIES
#ifdef OPENMP
/**
 * Split disconnected communities using Label Propagation Algorithm (LPA).
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vaff whether each vertex is affected, if pruning is enabled (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <bool PRUNE=false, class B, class G, class K>
inline void splitDisconnectedCommunitiesLpaOmpW(vector<K>& vcom, vector<B>& vaff, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    if constexpr (PRUNE) vaff[u] = B(1);
  }
  // Perform label propagation within each community.
  while (true) {
    size_t ndel = 0;
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:ndel)
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      if (PRUNE && !vaff[u]) continue;
      K d = vdom[u];
      K c = vcom[u];
      // Find the minimum label of all neighbors in the same community.
      x.forEachEdgeKey(u, [&](auto v) {
        if (vdom[v]==d) c = min(c, vcom[v]);
      });
      if (c==vcom[u]) continue;
      // Update the label of this vertex.
      vcom[u] = c;
      ++ndel;
      if constexpr (!PRUNE) continue;
      vaff[u] = B();
      x.forEachEdgeKey(u, [&](auto v) { if (vdom[v]==d && !vaff[v]) vaff[v] = B(1); });
    }
    if (ndel==0) break;
  }
}


/**
 * Split disconnected communities using DFS.
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vis vertex visited flags (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <class B, class G, class K>
inline void splitDisconnectedCommunitiesDfsOmpW(vector<K>& vcom, vector<B>& vis, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    vis[u]  = B();
  }
  // Perform DFS from an untouched vertex, within each community (each thread picks a community atomically).
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    int T = omp_get_num_threads();
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K d = vdom[u];
      K c = vcom[u];
      if (!belongsOmp(d, t, T) || vis[u]) continue;
      auto ft = [&](auto v) { return vdom[v]==d; };
      auto fp = [&](auto v) { vcom[v] = c; };
      dfsVisitedForEachU(vis, x, u, ft, fp);
    }
  }
}


/**
 * Split disconnected communities using BFS.
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vis vertex visited flags (scratch)
 * @param us per-thread start vertices for BFS (scratch)
 * @param vs per-thread frontier vertices for BFS (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <class B, class G, class K>
inline void splitDisconnectedCommunitiesBfsOmpW(vector<K>& vcom, vector<B>& vis, vector<vector<K>*>& us, vector<vector<K>*>& vs, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    vis[u]  = B();
  }
  // Perform DFS from an untouched vertex, within each community (each thread picks a community atomically).
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    int T = omp_get_num_threads();
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K d = vdom[u];
      K c = vcom[u];
      if (!belongsOmp(d, t, T) || vis[u]) continue;
      auto ft = [&](auto v, auto _) { return vdom[v]==d; };
      auto fp = [&](auto v, auto _) { vcom[v] = c; };
      (*us[t]).clear(); (*vs[t]).clear(); (*us[t]).push_back(u);
      bfsVisitedForEachU(vis, *us[t], *vs[t], x, ft, fp);
    }
  }
}


/**
 * Split disconnected communities using BFS, vertex-parallel approach.
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param cvis community visited flags (scratch)
 * @param vis vertex visited flags (scratch)
 * @param us start vertices for BFS (scratch)
 * @param vs frontier vertices for BFS (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <class B, class G, class K>
inline void splitDisconnectedCommunitiesBfsOmpW(vector<K>& vcom, vector<B>& cvis, vector<B>& vis, vector<B>& us, vector<B>& vs, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  size_t N = x.order();
  size_t NVIS = 0;
  B CDONE = B(1);
  // Clear visited flags.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    cvis[u] = B();
    vis[u]  = B();
    us[u]   = B();
    vs[u]   = B();
  }
  // Initiate BFS one vertex per community.
  #pragma omp parallel for schedule(static, 2048) reduction(+:NVIS)
  for (K u=0; u<S; ++u) {
    K d = vdom[u];
    if (cvis[d]) continue;
    cvis[d] = CDONE;
    vcom[u] = u;
    vis[u]  = B(1);
    us[u]   = B(1);
    ++NVIS;
  }
  // Perform BFS within each community to split disconnected communities.
  while (countValueOmp(vis, B(1)) < N) {
    // Perform BFS within each community, until no more vertices are visited.
    while (1) {
      size_t NVISDEL = 0;
      #pragma omp parallel for schedule(dynamic, 2048) reduction(+:NVISDEL)
      for (K u=0; u<S; ++u) {
        if (!us[u]) continue;
        K d = vdom[u];
        K c = vcom[u];
        x.forEachEdgeKey(u, [&](auto v) {
          if (vis[v] || vdom[v]!=d) return;
          vcom[v] = c;
          vis[v]  = B(1);
          vs[v]   = B(1);
          ++NVISDEL;
        });
      }
      swap(us, vs);
      fillValueOmpU(vs, B());
      NVIS += NVISDEL;
      if (NVISDEL==0) break;
    }
    // if (countValueOmp(vis, B(1))==N) break;
    // Reinitiate BFS from untouched vertices in each community.
    ++CDONE;
    #pragma omp parallel for schedule(static, 2048) reduction(+:NVIS)
    for (K u=0; u<S; ++u) {
      K d = vdom[u];
      if (vis[u] || cvis[d]==CDONE) continue;
      cvis[d] = CDONE;
      vcom[u] = u;
      vis[u]  = B(1);
      us[u]   = B(1);
      ++NVIS;
    }
  }
}
#endif
#pragma endregion
