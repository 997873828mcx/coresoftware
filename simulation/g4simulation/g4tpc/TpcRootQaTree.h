#ifndef G4TPC_TPCROOTQATREE_H
#define G4TPC_TPCROOTQATREE_H

#include <TTree.h>

namespace TpcRootQaTree
{
  // These QA trees are written to their selected output files only at End.
  // Keep every basket in memory until then: the current ROOT directory at
  // InitRun may belong to a different DST file. Writing just a tree header
  // elsewhere would leave its already-flushed baskets in that other file.
  // The caller owns the returned tree.
  inline TTree *make(const char *name, const char *title)
  {
    auto *tree = new TTree(name, title);
    tree->SetDirectory(nullptr);
    return tree;
  }

  inline void reset(TTree *&tree, const char *name, const char *title)
  {
    auto *replacement = make(name, title);
    delete tree;
    tree = replacement;
  }
}

#endif
