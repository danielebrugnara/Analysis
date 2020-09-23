#pragma once

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
#include "Hit.h"

enum LineType {
	EVENT_START,
	BLANK_EVENT,
	GAMMA_EVENT,
	POSITRON_EVENT,
	ELECTRON_EVENT,
	ION_EVENT,
	SPECTROMETER_EVENT,
	DIAMANT_HIT,
	AGATA_HIT,
	LABR_HIT,
	TARGET_HIT,
	HEADER_END,
	UNKNOWN_EVENT
};

enum LineType Classify(const std::string&, bool);
double Resolution_AGATA(const double&);


int readsimu(const std::string& inputfile, const std::string& rootoutput, bool SMEARING=false, int it_max=-1, bool gun=false)
{

	Particle* tmp_particle=0;
	std::vector <Particle> particles;
	std::vector <Hit> hits;

	auto *output = new TFile(rootoutput.c_str(),"recreate");
	auto *theTree = new TTree("SimulatedAgata","SimulatedAgata");
	TBranch *branch_particle = theTree->Branch("Particle", &particles);
	TBranch *branch_hit = theTree->Branch("Hits", &hits);
	std::string aLine;
	bool Line_skip=true;
    std::istringstream();

	auto *randm=new TRandom();
	double Resolution_DIAMANT=0.021;
	//double Resolution_AGATA=0.002;
	//double Resolution_AGATA=0.003/2.355;//numerator is FWHM

	double Energy, posX, posY, posZ;
	int nr, type, segment, crystal;

	int iteration =0;
	std::cout<<"\n";
	std::ifstream input ; 
	input.open(Form("%s",inputfile.c_str()));
	if(input.is_open()){
		while(input.good()){
			if(iteration>it_max && it_max !=-1) break;
			if(iteration != 0 && iteration % 10000==0) std::cout<<"\r" <<iteration;
			std::getline(input,aLine);
			if (Classify(aLine, gun)==HEADER_END) Line_skip=false;
			if(aLine.empty()||Line_skip) continue;
			std::istringstream oneline(aLine);
			switch(Classify(aLine, gun)){
				case EVENT_START:{//Start of event
							 theTree->Fill();
							 iteration++;
							 particles.clear();
							 hits.clear();
							 break;
						 }
				case BLANK_EVENT:{//Line whithout information
							 break;
						 }
				case GAMMA_EVENT:{//Gamma original event
							 oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							 if (tmp_particle!=nullptr){
								 particles.push_back(*tmp_particle);
								 delete tmp_particle;
								 tmp_particle = nullptr;
							 }
							 tmp_particle= new Particle(1, Energy, posX, posY, posZ, nr);                  
							 break;
						 }
				case POSITRON_EVENT:{//Positron original event
							    oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							    if (tmp_particle!= nullptr){
								    particles.push_back(*tmp_particle);
								    delete tmp_particle;
								    tmp_particle = nullptr;
							    }
							    tmp_particle= new Particle(98, Energy, posX, posY, posZ, nr);                  
							    break;
						    }
				case ELECTRON_EVENT:{//Electron original event
							    oneline >> type >> Energy >> posX >> posY >> posZ >> nr;
							    if (tmp_particle!= nullptr){
								    particles.push_back(*tmp_particle);
								    delete tmp_particle;
								    tmp_particle = nullptr;
							    }
							    tmp_particle= new Particle(97, Energy, posX, posY, posZ, nr);                  
							    break;
						    }
				case ION_EVENT:{//Ion original event
						       break;
					       }
				case SPECTROMETER_EVENT:{//Spectrometer original event
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
						       oneline >> crystal >> Energy >> posX >> posY >> posZ >> segment;
						       if (SMEARING && Energy>0){
                                   // std::cout << "before : " << Energy << std::endl;
							       Energy =randm->Gaus(Energy, Resolution_AGATA(Energy));
                                   // std::cout << "after : " << Energy << std::endl;
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
				case TARGET_HIT:{//Hit on cryo target
						      oneline >> crystal >> Energy >> posX >> posY >> posZ >> segment;
						      Hit tmp_hit(crystal, Energy, posX, posY, posZ, segment, 3);
						      hits.push_back(tmp_hit);
						      break;
					      }
				case HEADER_END:{//Line whithout information
							//  Line_skip=false;
							break;
						}
				case UNKNOWN_EVENT:{
							throw std::runtime_error("Something wrong in reading : "+aLine+"\n");
						}
				default:{}
			}            
		}
	}else{
	    std::cerr<<"unable to open\n";
	    return 0;
	};

	theTree->Write();
	output->Close();

    if (iteration == 0){
        return readsimu( inputfile, rootoutput, SMEARING, it_max, !gun);
    }
    return iteration-1;
}

enum LineType Classify(const std::string& Line, bool gun){
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
	if (!Line.compare(0, 5, "   -8")){
		if(!gun)
			return ION_EVENT;
		else
			return EVENT_START;
	}
	if(!Line.compare(0, 5, " -104"))
		return SPECTROMETER_EVENT;
	if (!Line.compare(0, 3, " 80"))
		return DIAMANT_HIT;
	if (!Line.compare(0, 5, "18002"))
		return LABR_HIT;
	if (!Line.compare(0, 5, " 4000"))
		return TARGET_HIT;
	if (!Line.compare(0, 1, "$"))
		return HEADER_END;
	if(([Line](){try{return (std::stoi(Line.substr(0, 5))>=0 && std::stoi(Line.substr(0, 5))<=100);}
		    catch(const std::invalid_argument & err){return false;}}()))
		return AGATA_HIT;
	else
		return UNKNOWN_EVENT;
}

double Resolution_AGATA(const double& energy){
//   1  p0           8.12318e-01           nan          -nan          -nan
//   2  p1           2.06151e-03           nan          -nan          -nan
    return (8.12318e-01+2.06151e-03*energy)/2.;
}
