#include<TFile.h>

#include<string>
#include<vector>
#include<map>

StripUptime(std::string file_name){
	TFile* file = new TFile(file_name.c_str());
	std::vector<std::string> mg = {"MG1","MG3","MG4","MG5","MG7","MG11"};
	std::vector<std::string> type = {"T","E"};
	std::vector<std::string> strip = {"X","Y"};
	std::map<std::string, std::map<std::string, std::vector<std::pair>>>> data;
	//for(const auto &it_type: type){
		for(const auto &it_mg: mg){
			for(const auto &it_strip: strip){
				std::string tmp_name = "pConf_SI_mStrip_"
							+it_type
							+"_"
							+it_mg
							+"_"
							+it_strip;
				TH2D* tmp_ptr = (TH2D*)file->Get(tmp_name);

				for( int ii=0; ii<tmp_ptr->GetXaxis()->GetNbins(); ++ii){
					double counts = 0;
					for( int jj=0; jj<tmp_ptr->GetYaxis()->GetNbins(); ++jj){
						int bin = tmp_ptr->GetBin(ii, jj);
						double cont = tmp_ptr->GetBinContent(bin);
						counts+=cont;
					}	
					data[it_type][it_mg][it_strip].push_back(counts);
				}

			}
		}
	//}
	for(const auto &it_type: type){
		for(const auto &it_mg: type){
			for(const auto &it_strip: strip){
				for( const auto & counts: data[it_type][it_mg][it_strip]){
							
				}
			}
		}
	}
}
