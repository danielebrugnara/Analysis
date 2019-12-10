#include "Main.h"

int main(int argc, char *argv[]){
    TApplication root_app("app",&argc,argv);
    //Analysis *my_analysis = nullptr;
    Analysis *my_analysis = nullptr;

    try{
        my_analysis = new Analysis();
        //my_analysis = new Analysis();

    }catch (const std::runtime_error& error){
        std::cerr << error.what() <<std::endl;
        return -1;
    }catch (const std::logic_error& error){
        std::cerr << "Caught error : " << error.what() << std::endl;
        return -1;

    }
    my_analysis->RunAnalysis(-1, 0);

    root_app.Run();
    return 0;
}