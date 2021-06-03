#include <iostream>
#include <string>

#include <sys/wait.h>
#include <unistd.h>
#include <dirent.h>

#include "RunSelector.h"

int main(int argc, char* argv[]){
    std::vector<std::string> files{
        //"./../../../DataAnalyzed/simu/46Ar3Hed47K_0keV_s12_0_0.root",
        //"./../../../DataAnalyzed/simu/46Ar3Hed47K_360keV_d32_0_0.root",
        //"./../../../DataAnalyzed/simu/46Ar3Hed47K_2020keV_f72_0_0.root",
    };
    DIR* dir;
    struct dirent* ent;
    std::string direct{"./../../../../DataAnalyzed/simu/"};
    if ((dir = opendir(direct.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            std::string tmpname(ent->d_name);
            std::string necessaryBegin = "46Ar3Hed47K_";
            std::string necessaryEnd = "um_analyzed.root" ;
            if(tmpname.compare(0, necessaryBegin.size(), necessaryBegin) != 0 ) continue;
            std::cout << tmpname << std::endl;
            if(not std::equal(necessaryEnd.rbegin(), necessaryEnd.rend(), tmpname.rbegin())) continue;
            printf("Adding: %s\n", ent->d_name);
            files.emplace_back(direct+tmpname);
        }
        closedir(dir);
    } else {
        /* could not open directory */
        perror("");
    }

    RunSelector run("../../../../DataAnalyzed/simu/sumsimu.root", "../../../../DataAnalyzed/sum.root", files);
    //RunSelector run("Data/anaout/anular.root", "../../../DataAnalyzed/sum.root");
    return 0;
}
