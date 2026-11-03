#include "TpcRootQaTree.h"

#include <TFile.h>
#include <TSystem.h>

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
  constexpr int entries = 30000;
  const std::array<const char *, 5> names = {
      "driftXY", "driftStepQA", "ionizationEventQA", "ionizationClusterQA", "ionizationSegmentQA"};

  void require(bool condition, const std::string &message)
  {
    if (!condition) throw std::runtime_error(message);
  }
}

// A small first event leaves all baskets in memory and can hide the bug.
// Exceed many basket boundaries while a DIFFERENT file is current, then write
// each tree to its final QA destination. No source file may be needed to read it.
int main()
{
  const std::string stem = "test-tpc-qa-" + std::to_string(gSystem->GetPid());
  const std::string foreign_name = stem + "-dst.root";
  const std::string qa_name = stem + "-qa.root";
  const std::string split_name = stem + "-split.root";
  try
  {
    {
      TFile foreign(foreign_name.c_str(), "RECREATE");
      require(!foreign.IsZombie(), "cannot create unrelated DST output");
      std::array<std::unique_ptr<TTree>, 5> trees;
      std::array<float, 5> values{};
      for (std::size_t i = 0; i < trees.size(); ++i)
      {
        foreign.cd();
        // Reset also supports a second InitRun without leaking the old tree.
        TTree *raw = TpcRootQaTree::make("previous_run", "previous run QA");
        TpcRootQaTree::reset(raw, names[i], "QA basket ownership regression");
        trees[i].reset(raw);
        require(trees[i]->GetDirectory() == nullptr, "QA tree attached to unrelated current file");
        auto *branch = trees[i]->Branch("value", &values[i], "value/F");
        require(branch != nullptr, "missing test branch");
        trees[i]->SetBasketSize("*", 1024);
      }
      for (int entry = 0; entry < entries; ++entry)
      {
        for (std::size_t i = 0; i < trees.size(); ++i)
        {
          values[i] = static_cast<float>(entry * 8 + static_cast<int>(i));
          require(trees[i]->Fill() > 0, "cannot fill QA tree");
        }
      }
      {
        TFile qa(qa_name.c_str(), "RECREATE");
        require(!qa.IsZombie(), "cannot create combined QA output");
        for (std::size_t i : {0UL, 2UL, 3UL, 4UL})
          require(trees[i]->Write() > 0, "cannot write combined QA tree");
      }
      {
        TFile split(split_name.c_str(), "RECREATE");
        require(!split.IsZombie(), "cannot create split step QA output");
        require(trees[1]->Write() > 0, "cannot write split step QA tree");
      }
      foreign.Close();
    }
    require(std::remove(foreign_name.c_str()) == 0, "cannot remove unrelated DST output");
    for (std::size_t i = 0; i < names.size(); ++i)
    {
      TFile file((i == 1 ? split_name : qa_name).c_str(), "READ");
      require(!file.IsZombie() && !file.TestBit(TFile::kRecovered), "QA output unreadable/recovered");
      auto *tree = file.Get<TTree>(names[i]);
      require(tree && tree->GetEntries() == entries, "QA rows missing");
      float value = 0.;
      require(tree->SetBranchAddress("value", &value) >= 0, "cannot bind QA branch");
      for (int entry = 0; entry < entries; ++entry)
      {
        require(tree->GetEntry(entry) > 0, "unreadable QA basket");
        require(value == static_cast<float>(entry * 8 + static_cast<int>(i)), "QA basket value changed");
      }
    }
    std::remove(qa_name.c_str());
    std::remove(split_name.c_str());
    std::cout << "PASS: 150000 QA rows survive basket boundaries and split outputs without the original DST\n";
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
