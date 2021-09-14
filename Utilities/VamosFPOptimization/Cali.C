#include <string>
#include <thread>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TProfile.h>
#include <TSpline.h>
#include <TLine.h>

#include <TROOT.h>
#include <TGClient.h>
#include <TSystem.h>

#include "./Minimizer.cxx"


std::string GetCommand();
std::string GetCommand2();
std::string GetCommand2h();
void SetAliases(TTree *);
std::pair<TH2D *, TH2D*>GenerateHisto(TTree *);
void DrawLines();
void LoadStates();
//void LoadProfiles();
double h= 0.01;
TSpline3* CaliTime();
double GetTime(double );

//std::map <std::string, TSpline3*> interpolations;
std::map <std::string, double> states;
std::vector<Minimizer*> minimizers;

const int bins = 60;

std::map <std::string, TProfile *> profiles1;
std::map <std::string, TProfile *> profiles2;
TSpline3* calibration = nullptr;
TH1D* shifts = nullptr;
TH1D* discrepancies1 = nullptr;
TH1D* discrepancies2 = nullptr;
int counter = 0;


struct Limits{
	double X_min;
	double X_max;
	double Y_min;
	double Y_max;
	Limits(double a, double b, double c, double d): 
			X_min(a), X_max(b), Y_min(c), Y_max(d){} 	
};

std::map <std::string, Limits*> limits;

void SetupLimits(){
	limits["14"] = new Limits(3.1955, 3.3535, -372, -140.8);
	limits["15"] = new Limits(3., 3.116, -204, -25.6);
	limits["16"] = new Limits(2.814, 2.937, -102.4, 125.6);
	limits["17"] = new Limits(2.659, 2.806, 33.6, 254.4);
	limits["18"] = new Limits(2.52, 2.69, 146, 359);
}

TProfile* GetProfile(TH2D* histo, double X_min, double X_max,
						double Y_min, double Y_max){	


	histo->GetXaxis()->SetRangeUser(X_min, X_max);
	histo->GetYaxis()->SetRangeUser(Y_min, Y_max);
	return histo->ProfileY(Form("profile_%s_%f_%f", histo->GetName(), X_min, Y_min));
}

void ComputeDiscrepancies()
{
	//delete discrepancies;
	discrepancies1 = new TH1D(Form("discr1_it%d", counter),Form("discr1_it%d", counter), bins, -500, 400);
	discrepancies2 = new TH1D(Form("discr2_it%d", counter),Form("discr2_it%d", counter), bins, -500, 400);
	++counter;
	for (int ii = 0; ii < bins; ++ii)
	{

		double point = discrepancies1->GetBinCenter(ii);
		double val = 0;
		double tmp;
		int cnt = 0;
		for (const auto &profile : profiles1)
		{
			//interpolation.second->Draw();
			if (point > profile.second->GetXaxis()->GetXmin() * 1.01 && point < profile.second->GetXaxis()->GetXmax() * 0.99)
			{
				tmp = profile.second->GetBinContent(profile.second->FindBin(point) ) - states[profile.first];
				//std::cout << "tmp: "<< tmp << " state["<<profile.first<<"] : " << states[profile.first] << std::endl;
				if (abs(tmp) > 0.5){
					//std::cout << "Skipping for large discrepance bin : " << ii <<std::endl;
					continue;
				}
				if (cnt == 0)
					val = tmp;
				else
					val += tmp;
				//std::cout << val << std::endl;
				++cnt;
			}
		}
		if (cnt == 0){
			//std::cout << "Warning, value zero for bin : " << ii <<std::endl;
			val = 0;
		}
		else
			val = val / cnt;
		discrepancies1->SetBinContent(ii, val);

		val = 0;
		cnt = 0;
		for (const auto &profile : profiles2)
		{
			//interpolation.second->Draw();
			if (point > profile.second->GetXaxis()->GetXmin() * 1.01 && point < profile.second->GetXaxis()->GetXmax() * 0.99)
			{
				tmp = profile.second->GetBinContent(profile.second->FindBin(point) ) - states[profile.first];
				//std::cout << "tmp: "<< tmp << " state["<<profile.first<<"] : " << states[profile.first] << std::endl;
				if (abs(tmp) > 0.5){
					//std::cout << "Skipping for large discrepance bin : " << ii <<std::endl;
					continue;
				}
				if (cnt == 0)
					val = tmp;
				else
					val += tmp;
				//std::cout << val << std::endl;
				++cnt;
			}
		}
		if (cnt == 0){
			//std::cout << "Warning, value zero for bin : " << ii <<std::endl;
			val = 0;
		}
		else
			val = val / cnt;
		discrepancies2->SetBinContent(ii, val);
	}
}


//int test = 911;

void InitializeMinimizers(){

	for (int ii=0; ii<bins;++ii){
		minimizers.push_back(new Minimizer(0., 5000, 0.01, 100, 0.99, h));
		//else minimizers.push_back(new Minimizer(discrepancies->GetBinContent(ii), true));
	}
}

void Minimize(){
	//delete shifts;
	shifts = new TH1D(Form("shifts_it%d", counter),Form("shifts_it%d", counter), bins, -500, 400);
	double value =0;
	for (int ii=0; ii<bins;++ii){
		value = minimizers.at(ii)->Step(std::make_pair(pow(discrepancies1->GetBinContent(ii), 2), pow(discrepancies2->GetBinContent(ii), 2)));
		if (abs(value)>20) std::cout << "SOME PROBLEM WITH BIN :" << ii <<"  value : " << value  <<std::endl;
		//std::cout << value << std::endl;
		shifts->SetBinContent(ii, value);
	}

	for (int ii=70; ii>-1; --ii){
		if (shifts->GetBinContent(ii)!=0) continue;
		shifts->SetBinContent(ii, shifts->GetBinContent(ii+1));
	}
	for (int ii=940; ii<bins; ++ii){
		if (shifts->GetBinContent(ii)!=0) continue;
		shifts->SetBinContent(ii, shifts->GetBinContent(ii-1));
	}

	for (int ii=1; ii<bins-1; ++ii){
        double current = shifts->GetBinContent(ii);
        double next = shifts->GetBinContent(ii+1);
        double prev = shifts->GetBinContent(ii-1);
        if (abs(current-0.5*(next+prev))>3 && ii>shifts->FindBin(-375) && ii<shifts->FindBin(350)){
            std::cout << "discrepancy out of line in bin " << ii << std::endl;
            shifts->SetBinContent(ii, 0.5*(next+prev)); 
        }
    
    }
	calibration = new TSpline3(shifts);
}

void Cali(const char* file_name){

	TFile * file = new TFile(file_name, "read");
	TTree * tree =(TTree*) file->Get("PhysicsTree");
	SetAliases(tree);
	SetupLimits();
	LoadStates();
	//std::vector<std::thread> threads;
	bool first = true;

	unsigned int iteration = 0;
	try{
		while (++iteration<19)
		{
			std::pair<TH2D *, TH2D*> histo = GenerateHisto(tree);


			for (const auto &limit : limits)
			{
				profiles1[limit.first] = GetProfile(	histo.first,
									limit.second->X_min, limit.second->X_max,
									limit.second->Y_min, limit.second->Y_max);
				histo.first->GetXaxis()->UnZoom();
				histo.first->GetYaxis()->UnZoom();

				profiles2[limit.first] = GetProfile(	histo.second,
									limit.second->X_min, limit.second->X_max,
									limit.second->Y_min, limit.second->Y_max);
				histo.second->GetXaxis()->UnZoom();
				histo.second->GetYaxis()->UnZoom();

			}
			ComputeDiscrepancies();

			if (first){ InitializeMinimizers(); first= false;};
			Minimize();
			//delete calibration;
			//shifts->Draw("");
			//while (gROOT->FindObject("cvvv") != NULL){
			//	gSystem->Sleep(200);
			//	gClient->HandleInput();
			//	gSystem->ProcessEvents();
			//}
		}
	}
	catch(const std::runtime_error err){
		std::cout << "finished iterating :" << err.what() << std::endl;

	} 
	//TCanvas *final_cv = new TCanvas("FINAL_CANVAS", "FINAL_CANVAS");
    //tree->Draw(Form("Xf:Mass_Q>>FINAL_HISTO(2000, 2.4, 3.4, %i, -400, 400)", bins), "IC[0]>0.1&&IC[1]>0.1", "histo col");

	//for (int ii=0; ii<bins;++ii){
	//std::cout << "bin : " << ii << "  shift : " << shifts->GetBinContent(ii) <<std::endl;
	//}

	//DrawLines();
	new TCanvas("finalcanvas", "finalcanvas");
	calibration->Draw();
}

void DrawLines(){
    TLine *line1 = new TLine(3.28, -500, 3.28, 400);
    TLine *line2 = new TLine(3.07, -500, 3.07, 400);
    TLine *line3 = new TLine(2.88, -500, 2.88, 400);
    TLine *line4 = new TLine(2.71, -500, 2.71, 400);
    TLine *line5 = new TLine(2.56, -500, 2.56, 400);

    line1->Draw();
    line2->Draw();
    line3->Draw();
    line4->Draw();
    line5->Draw();

}

std::pair<TH2D *, TH2D*>GenerateHisto(TTree *tree){
    TCanvas* cv = new TCanvas(Form("PROGRESS_%i", counter));
    tree->Draw(Form("Xf:Mass_Q>>hh1_%i(2000, 2.4, 3.4, %i, -500, 400)", counter, bins), "IC[0]>0.1&&IC[1]>0.1", "histo col");
    tree->Draw(Form("Xf:Mass_Qh>>hh2_%i(2000, 2.4, 3.4, %i, -500, 400)", counter, bins), "IC[0]>0.1&&IC[1]>0.1", "histo col");
	TH2D *hh1 = (TH2D *)gDirectory->Get(Form("hh1_%i", counter));
	DrawLines();
	TH2D *hh2 = (TH2D *)gDirectory->Get(Form("hh2_%i", counter));
	DrawLines();
	return std::make_pair(hh1, hh2);
}

std::string GetCommand(){
	std::string cmd = "540.5*(AGAVA_VAMOSTS<104753375647998)+537.9*(AGAVA_VAMOSTS>=104753375647998) -2.*T_FPMW_CATS2_C";
	std::ifstream cali_file("FP_Time.cal");
	std::string line;
	std::string min;
	std::string max;
	std::string shift;
	if (!cali_file) std::cerr <<"Error Opening cali file\n";
	while (std::getline(cali_file, line)){
		std::stringstream str;
		str << line;
		str >> min >> max >> shift;
		cmd.append("+"+shift+"*"+"(Xf>"+min+")*(Xf<"+max+")");
	}
	//cmd.append("+GetTime(Xf)");

	std::cout << cmd <<std::endl;
	return cmd;
}

//std::string GetCommand2(){
//	std::string cmd = "540.5*(AGAVA_VAMOSTS<104753375647998)+537.9*(AGAVA_VAMOSTS>=104753375647998) -2.*T_FPMW_CATS2_C";
//	std::ifstream cali_file("FP_Time.cal");
//	std::string line;
//	std::string min;
//	std::string max;
//	std::string shift;
//	if (!cali_file) std::cerr <<"Error Opening cali file\n";
//	while (std::getline(cali_file, line)){
//		std::stringstream str;
//		str << line;
//		str >> min >> max >> shift;
//		cmd.append("+("+shift+"+"+std::to_string(h)+")*"+"(Xf>"+min+")*(Xf<"+max+")");
//	}
//	cmd.append("+GetTime(Xf)");
//	return cmd;
//}
std::string GetCommand2(){
	std::string cmd = "Path/cos(Pf/1000.)+(760.-752.81)/(cos(Pf/1000.)*cos(Tf/1000.))";
	cmd.append("+GetTime(Xf)");
	return cmd;
}
std::string GetCommand2h(){
	std::string cmd = "Path/cos(Pf/1000.)+(760.-752.81)/(cos(Pf/1000.)*cos(Tf/1000.))";
	cmd.append("+GetTime(Xf)+"+std::to_string(h));
	return cmd;
}


void SetAliases(TTree *tree){
	tree->SetAlias("IC0", "IC[0]*1.04");
	tree->SetAlias("IC1", "IC[1]*1.03");
	tree->SetAlias("IC2", "IC[2]*0.98");
	tree->SetAlias("IC3", "IC[3]*1.");
	tree->SetAlias("IC4", "IC[4]*1.");
	tree->SetAlias("IC5", "IC[5]*1.");
	tree->SetAlias("IC6", "IC[6]*1.");
	tree->SetAlias("En", "IC0+IC1+IC2+IC3+IC4+IC5+IC6");
	tree->SetAlias("dEn0", "IC0");
	tree->SetAlias("dEn1", "IC1");
	std::string cmd = GetCommand(); 
	std::string cmd2 = GetCommand2(); 
	std::string cmd2h = GetCommand2h(); 
	//std::string cmd2 = GetCommand2(); 
   	tree->SetAlias("Time",cmd.c_str());	
	tree->SetAlias("Dis", cmd2.c_str());
	tree->SetAlias("Dish", cmd2h.c_str());
   	//tree->SetAlias("Timeh",cmd2.c_str());	
	tree->SetAlias("Vel", "Dis/Time");
	tree->SetAlias("Velh", "Dish/Time");
	tree->SetAlias("Bet", "Vel/29.9792");
	tree->SetAlias("Beth", "Velh/29.9792");
	tree->SetAlias("Gamm", "1./sqrt(1-Bet*Bet)");
	tree->SetAlias("Gammh", "1./sqrt(1-Beth*Beth)");
	tree->SetAlias("Mass_Q", "Brho/3.107/Bet/Gamm");
	tree->SetAlias("Mass_Qh", "Brho/3.107/Beth/Gammh");
}

void LoadStates(){
	states["18"] = 46./18. ;//2.56;
	states["17"] = 46./17. ;//2.71;
	states["16"] = 46./16. ;//2.88;
	states["15"] = 46./15. ;//3.07;
	states["14"] = 46./14. ;//3.28;
}

double GetTime(double Xf){
	if (calibration){ 
		return calibration->Eval(Xf);
	}
	return 0;
}
