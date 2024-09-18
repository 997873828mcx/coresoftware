// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_BCORANGECHECK_H
#define FFARAWMODULES_BCORANGECHECK_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class BcoRangeCheck : public SubsysReco
{
 public:
  BcoRangeCheck(const std::string &name = "BcoRangeCheck");

  ~BcoRangeCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;
  //  int ResetEvent(PHCompositeNode *topNode) override;

  void MyEvtNode(const std::string &name) { m_EvtNodeName = name; }

 private:
  std::string m_EvtNodeName = "INTTRAWHIT";
  std::set<uint64_t> bclk_seen;
  std::map<uint64_t, int> diffcnt;
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H