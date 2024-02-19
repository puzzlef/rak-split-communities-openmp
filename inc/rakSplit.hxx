#pragma once
#include <utility>
#include <tuple>
#include <vector>
#include <cstdint>
#include "_main.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::move;
using std::get;




#pragma region METHODS
#pragma region ENVIRONMENT SETUP
#ifdef OPENMP
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param o rak options
 * @param fi initialzing community membership (vcom)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom)
 * @param fa is vertex allowed to be updated? (u)
 * @returns rak result
 */
template <int SPLIT=0, bool DYNAMIC=false, class FLAG=char, class G, class FI, class FM, class FA>
inline auto rakSplitLastInvokeOmp(const G& x, const RakOptions& o, FI fi, FM fm, FA fa) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using F = FLAG;
  int l = 0;
  int T = omp_get_max_threads();
  // Get graph properties.
  size_t S = x.span();
  size_t N = x.order();
  // Allocate buffers.
  vector<F> vaff(S);  // Affected vertex flag
  vector<K> vcom;     // Community membership
  vector<vector<K>*> vcs(T);    // Hashtable keys
  vector<vector<W>*> vcout(T);  // Hashtable values
  if (!DYNAMIC) vcom.resize(S);
  rakAllocateHashtablesW(vcs, vcout, S);
  // Data structures for splitting disconnected communities.
  vector<vector<K>*> us(T), vs(T);        // Per-thread start, frontier vertices for BFS
  if (SPLIT) {
    for (int t=0; t<T; ++t) {
      us[t] = new vector<K>();
      vs[t] = new vector<K>();
      (*us[t]).reserve(4*S/T);
      (*vs[t]).reserve(4*S/T);
    }
  }
  // Perform RAK algorithm.
  float tm = 0, ti = 0, ts = 0;  // Time spent in different phases
  float t  = measureDuration([&]() {
    // Initialize community membership.
    ti += measureDuration([&]() { fi(vcom); });
    // Mark affected vertices.
    tm += measureDuration([&]() { fm(vaff, vcs, vcout, vcom); });
    // Perform iterations.
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationOmpW(vcom, vaff, vcs, vcout, x, fa); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
    ts += measureDuration([&]() {
      if (SPLIT==1)      { splitDisconnectedCommunitiesLpaOmpW<false>(vcom, vaff, x, ucom);  swap(ucom, vcom); }
      else if (SPLIT==2) { splitDisconnectedCommunitiesLpaOmpW<true> (vcom, vaff, x, ucom);  swap(ucom, vcom); }
      else if (SPLIT==3) { splitDisconnectedCommunitiesDfsOmpW(vcom, vaff, x, ucom);         swap(ucom, vcom); }
      else if (SPLIT==4) { splitDisconnectedCommunitiesBfsOmpW(vcom, vaff, us, vs, x, ucom); swap(ucom, vcom); }
    });
  }, o.repeat);
  if (SPLIT) {
    for (int t=0; t<T; ++t) {
      delete us[t];
      delete vs[t];
    }
  }
  rakFreeHashtablesW(vcs, vcout);
  return RakResult<K>(vcom, l, t, tm/o.repeat, ti/o.repeat, ts/o.repeat, countValueOmp(vaff, F(1)));
}
#endif
#pragma endregion




#pragma region STATIC
#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static RAK.
 * @param x original graph
 * @param o rak options
 * @returns rak result
 */
template <int SPLIT=0, class FLAG=char, class G>
inline auto rakSplitLastStaticOmp(const G& x, const RakOptions& o={}) {
  auto fi = [&](auto& vcom) { rakInitializeOmpW(vcom, x); };
  auto fm = [ ](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) { fillValueOmpU(vaff, FLAG(1)); };
  auto fa = [ ](auto u) { return true; };
  return rakSplitLastInvokeOmp<SPLIT, false, FLAG>(x, o, fi, fm, fa);
}
#endif
#pragma endregion
#pragma endregion
