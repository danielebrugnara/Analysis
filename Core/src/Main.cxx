#include "Main.h"

int main(int argc, char* argv[]) {
    TApplication root_app("app", &argc, argv);
    Analysis* my_analysis = nullptr;

    int n_threads {1};
    if (argc == 2) {
        n_threads = std::stoi(argv[1]);
        std::cout << "Number of threads : " << n_threads << std::endl;
    }
    try {
        std::cout << "Starting Analysis\n";
        my_analysis = new Analysis(n_threads);
        my_analysis->RunAnalysis();
    } catch (const std::runtime_error& error) {
        std::cerr << "Caught error : " << error.what() << std::endl;
        return 3;
    } catch (const std::logic_error& error) {
        std::cerr << "Caught error : " << error.what() << std::endl;
        return 2;
    } catch (int error){
        std::cerr << "Caught error : " << error << std::endl; 
        return 1;
    }
    std::cout << "Analysis completed succesfully!\n";
    return 0;
}
