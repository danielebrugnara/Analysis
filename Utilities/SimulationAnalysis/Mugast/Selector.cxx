#define Selector_cxx
// The class definition in Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//


#include "Selector.h"
#include <TH2.h>
#include <TStyle.h>

void Selector::Begin(TTree * /*tree*/)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();

	angdistr = new TH1D("angdistr", "Angular distribution", 1000, 0, 3.1415);
	fOutput->Add(angdistr);
	ex = new TH1D("ex", "ex", 1000, -10, 10);
	fOutput->Add(ex);
	ex_theta = new TH2D("ex_theta", "ex vs theta", 1000, 0, 3.1415, 1000, -10, 10);
	fOutput->Add(ex_theta);
	theta_cm = new TH1D("theta_cm", "theta_cm", 1000, 0, 3.1415);
	fOutput->Add(theta_cm);
	kinematicline = new TH2D("kinematicline", "Kinematic line", 1000, 0, 3.1415, 1000, -20, 20);
	fOutput->Add(kinematicline);

	for (const auto & it: MG_nr){
		strip_E[it] = new TH2D(Form("MG%i_strip_E", it), Form("MG%i_strip_E", it), 129, 0, 129, 1000, 0, 30);
	}

}

Bool_t Selector::Process(Long64_t entry)
{
	// The Process() function is called for each entry in the tree (or possibly
	// keyed object in the case of PROOF) to be processed. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	// When processing keyed objects with PROOF, the object is already loaded
	// and is available via the fObject pointer.
	//
	// This function should contain the \"body\" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	//
	// Use fStatus to set the return value of TTree::Process().
	//
	// The return value is currently not used.

	fReader.SetLocalEntry(entry);
	//std::cout << *Nev << std::endl;
	//
	//
	//if (Ex[0]>-100)
	//	std::cout << Ex[0] << std::endl;
	for (int ii=0; ii<20;++ii){
		if (X[ii]>-100){
			position.SetXYZ(X[ii],Y[ii], Z[ii]);
			angdistr->Fill(position.Theta());
			ex->Fill(Ex[ii]);
			ex_theta->Fill(position.Theta(),Ex[ii]);
            theta_cm->Fill(ThetaCM[ii]);
			kinematicline->Fill(position.Theta(),ELab[ii]);
		}
	}


	return kTRUE;
}

void Selector::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
	TFile *top = new TFile("output.root", "recreate");
	TIter iter(fOutput);
	TObject *obj;
	while ((obj = iter()))
	{
		obj->Write();
	}
	top->Close();

}

void Selector::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.

}
