#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif




template <class G, class K, class V>
double getModularity(const G& x, const RakResult<K>& a, V M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, V(1));
}


template <class G>
void runExperiment(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  auto M = edgeWeight(x)/2;
  auto Q = modularity(x, M, 1.0f);
  printf("[%01.6f modularity] noop\n", Q);
  RakOptions o = {repeat};

  for (int i=0, f=10; f<=10000; f*=i&1? 5:2, ++i) {
    float tolerance = 1.0f / f;
    // Find RAK using a single thread (non-strict).
    auto ak = rakSeqStatic<false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] rakSeqStatic       {tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), tolerance);
    // Find RAK using a single thread (strict).
    auto al = rakSeqStatic<true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] rakSeqStaticStrict {tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), tolerance);
  }
}


int main(int argc, char **argv) {
  using K = int;
  using V = TYPE;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x;  // V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW<true>(x, file); println(x);
  auto y = symmetricize(x); print(y); printf(" (symmetricize)\n");
  // auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  runExperiment(y, repeat);
  printf("\n");
  return 0;
}
