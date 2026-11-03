//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EvtGenExtDecayerPhysics.cc, 2022/08/09 23:49:39 mccumber Exp $
//
/// \file derived from eventgenerator/pythia/decayer6/src/P6DExtDecayerPhysics.cc
/// \brief Implementation of the EvtGenExtDecayerPhysics class
///
/// \author Zhaozhong Shi; LANL, Los Alamos

#include "EvtGenExtDecayerPhysics.hh"

#include <Geant4/G4Decay.hh>
#include <Geant4/G4DecayTable.hh>
#include <Geant4/G4Exception.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTable.hh>  // for G4ParticleTable::G4PTblDi...
#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4ProcessVector.hh>  // for G4ProcessVector
#include <Geant4/G4String.hh>         // for G4String
#include <Geant4/G4VPhysicsConstructor.hh>
#include <Geant4/G4VProcess.hh>  // for G4VProcess
#include <Geant4/G4Version.hh>

#include <boost/io/ios_state.hpp>

#include <cstddef>  // for size_t
#include <fstream>
#include <iostream>  // for operator<<, endl, basic_o...
#include <set>
#include <sstream>
#include <string>  // for operator<<

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// hack for Geant 10.03.p03
// Geant 10.02.p02 has this defined
#ifndef aParticleIterator
#define aParticleIterator ((subInstanceManager.offset[g4vpcInstanceID])._aParticleIterator)
#endif

EvtGenExtDecayerPhysics::EvtGenExtDecayerPhysics(const G4String& name)
  : G4VPhysicsConstructor(name)
{
  /// Standard constructor
}

//
// protected methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EvtGenExtDecayerPhysics::ConstructProcess()
{
  // EvtGen model initialization changed the behavior of std::cout float print. Recover the default float print format
  boost::io::ios_all_saver ias_cout(std::cout);
  boost::io::ios_all_saver ias_cerr(std::cerr);

  /// Loop over all particles instantiated and add external decayer
  /// to all decay processes if External decayer is set

  // Create Geant4 external decayer
  extDecayer = new G4EvtGenDecayer();
  extDecayer->SetVerboseLevel(0);

  if (!m_decayFile.empty())
  {
    std::ifstream f_local(m_decayFile);

    if (f_local.good())
    {
      extDecayer->SetDecayTable(m_decayFile, false);
    }
    else
    {
      const char* CALIBRATIONROOT = getenv("CALIBRATIONROOT");
      if (!CALIBRATIONROOT)
      {
        exit(1);
      }
      std::string DecayFileCalibrationRepo = std::string(CALIBRATIONROOT) + "/EvtGen/" + m_decayFile;

      std::ifstream f_calib(DecayFileCalibrationRepo);

      if (f_calib.good())
      {
        extDecayer->SetDecayTable(DecayFileCalibrationRepo, false);
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << " : Fatal error. Decay file can not be found at "
                  << m_decayFile << " nor at " << DecayFileCalibrationRepo << std::endl;

        exit(1);
      }
    }
  }

  for (const int pdgCode : m_suppressPhotosPdgCodes)
  {
    extDecayer->SuppressPhotosForPdg(pdgCode);
    std::cout << "Suppressing PHOTOS radiation for EvtGen mother PDG "
              << pdgCode << std::endl;
  }

  aParticleIterator->reset();
  int decayer_used = 0;
  std::set<int> unresolvedOverridePdgCodes = m_overrideNativeDecayPdgCodes;
  while ((*aParticleIterator)())
  {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    const int pdgCode = particle->GetPDGEncoding();
    const bool overrideNativeDecay =
        m_overrideNativeDecayPdgCodes.find(pdgCode) !=
        m_overrideNativeDecayPdgCodes.end();
    bool foundDecayProcess = false;

    if (verboseLevel > 1)
    {
      std::cout << "Setting ext decayer for: "
                << aParticleIterator->value()->GetParticleName()
                << std::endl;
    }

    G4ProcessVector* processVector = pmanager->GetProcessList();
#if G4VERSION_NUMBER >= 1060
    for (size_t i = 0; i < processVector->length(); i++)
#else
    for (G4int i = 0; i < processVector->length(); i++)
#endif
    {
      G4Decay* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
      if (decay)
      {
        foundDecayProcess = true;
        // The extDecayer will be deleted in G4Decay destructor
        // increment counter in case we want to print out stats
        // for whatever reason (non null means it is used and
        // must not be deleted)
        decay->SetExtDecayer(extDecayer);
        decayer_used++;
      }
    }

    if (overrideNativeDecay && foundDecayProcess)
    {
      // G4Decay invokes its external decayer only when GetDecayTable() is
      // null. Delete and clear this explicitly selected particle's native
      // table; all other particle definitions are left unchanged.
      G4DecayTable* nativeDecayTable = particle->GetDecayTable();
      particle->SetDecayTable(nullptr);
      delete nativeDecayTable;
      unresolvedOverridePdgCodes.erase(pdgCode);

      std::cout << "Using EvtGen instead of the native decay table for "
                << particle->GetParticleName()
                << " (PDG " << pdgCode << ")" << std::endl;
    }
  }

  if (!unresolvedOverridePdgCodes.empty())
  {
    std::ostringstream message;
    message << "Could not override the native decay table with EvtGen for PDG code(s):";
    for (const int pdgCode : unresolvedOverridePdgCodes)
    {
      message << ' ' << pdgCode;
    }
    message << ". The particle definition or its G4Decay process is missing.";
    G4Exception("EvtGenExtDecayerPhysics::ConstructProcess",
                "EvtGenExtDecayerPhysics001", FatalException,
                message.str().c_str());
  }

  // If the extDecayer isn't used for this particle we need to delete it here
  // as far as I see this never happens but this makes coverity happy
  if (!decayer_used)
  {
    delete extDecayer;
  }

  if (verboseLevel > 0)
  {
    std::cout << "External decayer physics constructed." << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
