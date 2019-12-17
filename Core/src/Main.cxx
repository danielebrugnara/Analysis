#include "Main.h"

int main(int argc, char *argv[]){
    TApplication root_app("app",&argc,argv);
    Analysis *my_analysis = nullptr;

    int n_threads = 1;
    if (argc==2) {
        n_threads = stoi(argv[1]);
        std::cout << "Number of threads : " << n_threads << std::endl;
    }
    try{
        my_analysis = new Analysis(n_threads);
        my_analysis->RunAnalysis();
        //root_app.Run();

    }catch (const std::runtime_error& error){
        std::cerr << error.what() <<std::endl;
        return -1;
    }catch (const std::logic_error& error){
        std::cerr << "Caught error : " << error.what() << std::endl;
        return -1;

    }
    return 0;
}
