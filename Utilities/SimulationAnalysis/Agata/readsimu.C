#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <map>
#include <vector>

#include <TRandom.h>
#include <TF1.h>

#include "Particle.h"
#include "loader.C"

enum LineType {
	EVENT_START,
	BLANK_EVENT,
	GAMMA_EVENT,
	POSITRON_EVENT,
	ELECTRON_EVENT,
	ION_EVENT,
	HEADER_END,
	DIAMANT_HIT,
	AGATA_HIT,
	LABR_HIT
};

enum LineType Classify(std::string Line);


void readsimu(std::string inputfile, std::string rootoutput, bool SMEARING=false, int it_max=-1)
{

	Particle* tmp_particle=0;
	std::vector <Particle> particles;
	std::vector <Hit> hits;

	TFile *output = new TFile(rootoutput.c_str(),"recreate");
	TTree *theTree = new TTree("SimulatedAgata","SimulatedAgata");
	TBranch *branch_particle = theTree->Branch("Particle", &particles);
	TBranch *branch_hit = theTree->Branch("Hits", &hits);
	std::string aLine;
	bool Line_skip=true;
	std::istringstream oneline;

	TRandom *randm=new TRandom();
	double Resolution_DIAMANT=0.021;
	double Resolution_AGATA=0.002;

	double Energy, posX, posY, posZ;
	int nr, type, segment, crystal;

	int iteration =0;
	std::cout<<"\n";
	std::ifstream input ; 
	input.open(Form("%s",inputfile.c_str()));
	if(input.is_open()){
		while(input.good()){
			if(iteration>it_max && it_max !=-1) break;
			if(iteration % 10000==0) std::cout<<"\r" <<iteration;
			std::getline(input,aLine);
			if (Classify(aLine)==HEADER_END) Line_skip=false;
			if(aLine.empty()||Line_skip) continue;
			std::istringstream oneline(aLine);
			switch(Classify(aLine)){
				case EVENT_START:{//Start of event
							 theTree->Fill();
							 iteration++;
							 particles.clear();
							 hits.clear();
							 break;
						 }
				case HEADER_END:{//Line whithout information
							//  Line_skip=false;
							break;
						}
				case BLANK_EVENT:{//Line whithout information
							 break;
						 }
				case GAMMA_EVENT:{//Gamma original event
							 oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							 if (tmp_particle!=0){
								 particles.push_back(*tmp_particle);
								 delete tmp_particle;
							 }
							 tmp_particle= new Particle(1, Energy, posX, posY, posZ, nr);                  
							 break;
						 }
				case ELECTRON_EVENT:{//Electron original event
							    oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							    if (tmp_particle!=0){
								    particles.push_back(*tmp_particle);
								    delete tmp_particle;
							    }
							    tmp_particle= new Particle(97, Energy, posX, posY, posZ, nr);                  
							    break;
						    }
				case POSITRON_EVENT:{//Positron original event
							    oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							    if (tmp_particle!=0){
								    particles.push_back(*tmp_particle);
								    delete tmp_particle;
							    }
							    tmp_particle= new Particle(98, Energy, posX, posY, posZ, nr);                  
							    break;
						    }
				case ION_EVENT:{//Ion original event
						       break;
					       }
				case DIAMANT_HIT:{//Hit on Diamant
							 oneline >> crystal >> Energy >> posX >> posY >> posZ >> segment;
							 if (SMEARING && Energy>0){
								 Energy =randm->Gaus(Energy, Resolution_DIAMANT*Energy);
							 }
							 Hit tmp_hit(crystal, Energy, posX, posY, posZ, segment, 1);
							 hits.push_back(tmp_hit);
							 break;
						 }
				case AGATA_HIT:{//Hit on Agata
						       oneline >> segment >> Energy >> posX >> posY >> posZ >> crystal;
						       if (SMEARING && Energy>0){
							       Energy =randm->Gaus(Energy, Resolution_AGATA*Energy);
						       }
						       Hit tmp_hit(crystal, Energy, posX, posY, posZ, segment, 0);
						       hits.push_back(tmp_hit);
						       break;
					       }
				case LABR_HIT:{//Hit on LABR
						      oneline >> crystal >> Energy >> posX >> posY >> posZ >> segment;
						      Hit tmp_hit(crystal, Energy, posX, posY, posZ, segment, 2);
						      hits.push_back(tmp_hit);
						      break;
					      }
			}            
		}
	}else{std::cout<<"unable to open\n";};

	theTree->Write();
	output->Close();
}

enum LineType Classify(std::string Line){
	if(!Line.compare(0, 4, "-100"))
		return EVENT_START;
	if(!Line.compare(0, 5, " -101")||!Line.compare(0, 5, " -102"))
		return BLANK_EVENT;
	if (!Line.compare(0, 5, "   -1"))
		return GAMMA_EVENT;
	if (!Line.compare(0, 5, "  -98"))
		return POSITRON_EVENT;
	if (!Line.compare(0, 5, "  -97"))
		return ELECTRON_EVENT;
	if (!Line.compare(0, 5, "  -8"))
		return ION_EVENT;
	if (!Line.compare(0, 5, "18002"))
		return LABR_HIT;
	if (!Line.compare(0, 1, "$"))
		return HEADER_END;
	if (!Line.compare(0, 3, " 80"))
		return DIAMANT_HIT;
	else
		return AGATA_HIT;
}






