/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <map>
#include <vector>

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include "TApplication.h"
#include "TROOT.h"

#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TDatabasePDG.h"

using namespace std;

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesG4me";
  stringstream message;
  TFile *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch eventStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  GenParticle *gen;
  HepMCEvent *element, *eve;
  Candidate *candidate;
  Int_t pdgCode;

  const Double_t c_light = 2.99792458E8;

  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Long64_t eventCounter, numberOfEvents;

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file"
         << " output_file"
         << " input_file(s)" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in G4me format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputFile = TFile::Open(argv[2], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't open " << argv[2] << endl;
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());

    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);

    TChain *chain = new TChain("Delphes");

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();
    
    static const int kMaxTracks = 32768; //1024;//1048576;
    Candidate *candidates[kMaxTracks];
    struct Tracks_t {
      int    n;
      char   proc[kMaxTracks];
      char   sproc[kMaxTracks];
      int    status[kMaxTracks];
      int    parent[kMaxTracks];
      int    pdg[kMaxTracks];
      double vt[kMaxTracks];
      double vx[kMaxTracks];
      double vy[kMaxTracks];
      double vz[kMaxTracks];
      double  e[kMaxTracks];
      double px[kMaxTracks];
      double py[kMaxTracks];
      double pz[kMaxTracks];
    } tracks;
    TTree *tree_tracks;
    auto pdgdb = TDatabasePDG::Instance();      
    
    for(i = 3; i < argc && !interrupted; ++i)
    {
      cout << "** Reading " << argv[i] << endl;

      inputFile = TFile::Open(argv[i]);

      if(inputFile == NULL)
      {
        message << "can't open " << argv[i] << endl;
        throw runtime_error(message.str());
      }

      tree_tracks = (TTree *)inputFile->Get("Tracks");
      tree_tracks->SetBranchAddress("n"      , &tracks.n);
      tree_tracks->SetBranchAddress("proc"   , &tracks.proc);
      tree_tracks->SetBranchAddress("sproc"  , &tracks.sproc);
      tree_tracks->SetBranchAddress("status" , &tracks.status);
      tree_tracks->SetBranchAddress("parent" , &tracks.parent);
      tree_tracks->SetBranchAddress("pdg"    , &tracks.pdg);
      tree_tracks->SetBranchAddress("vt"     , &tracks.vt);
      tree_tracks->SetBranchAddress("vx"     , &tracks.vx);
      tree_tracks->SetBranchAddress("vy"     , &tracks.vy);
      tree_tracks->SetBranchAddress("vz"     , &tracks.vz);
      tree_tracks->SetBranchAddress("e"      , &tracks.e);
      tree_tracks->SetBranchAddress("px"     , &tracks.px);
      tree_tracks->SetBranchAddress("py"     , &tracks.py);
      tree_tracks->SetBranchAddress("pz"     , &tracks.pz);
      numberOfEvents = tree_tracks->GetEntries();

      if(numberOfEvents <= 0) continue;

      // ExRootProgressBar progressBar(numberOfEvents - 1);
      ExRootProgressBar progressBar(-1);

      // Loop over all objects
      eventCounter = 0;
      modularDelphes->Clear();
      treeWriter->Clear();
      for(Int_t entry = 0; entry < numberOfEvents && !interrupted; ++entry)
      {

	tree_tracks->GetEntry(entry);

        for(Int_t j = 0; j < tracks.n; j++)
        {

	  auto ppdg = pdgdb->GetParticle(tracks.pdg[j]);
	  if (!ppdg) continue;
	  
          candidate = factory->NewCandidate();

	  // mm/c 
	  
          candidate->Momentum.SetPxPyPzE(tracks.px[j], // [GeV]
					 tracks.py[j], // [GeV]
					 tracks.pz[j], // [GeV]
					 tracks.e[j]); // [GeV]
	  
          candidate->Position.SetXYZT(tracks.vx[j] * 10.,              // [cm] -> [mm]
				      tracks.vy[j] * 10.,              // [cm] -> [mm]
				      tracks.vz[j] * 10.,              // [cm] -> [mm]
				      tracks.vt[j] * 1.e-6 * c_light); // [ns] -> [mm/c]

          candidate->PID = tracks.pdg[j];
          candidate->Status = 1;

          candidate->M1 = tracks.parent[j];
          candidate->M2 = -1;

          candidate->D1 = -1;
          candidate->D2 = -1;

          candidate->Charge = ppdg->Charge() / 3.;
          candidate->Mass = ppdg->Mass();

	  allParticleOutputArray->Add(candidate);
	  stableParticleOutputArray->Add(candidate);
	  candidates[j] = candidate;

	  // if has parent, set parent daughter label
	  auto imoth = candidate->M1;
	  if (imoth == -1) continue;
	  if (candidates[imoth]->D1 == -1) {
	    candidates[imoth]->D1 = j;
	    continue;
	  }
	  candidates[imoth]->D2 = j;
	  
	}

        modularDelphes->ProcessTask();

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();

        progressBar.Update(eventCounter, eventCounter);
        ++eventCounter;
      }

      progressBar.Update(eventCounter, eventCounter, kTRUE);
      progressBar.Finish();

      inputFile->Close();

    }

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
