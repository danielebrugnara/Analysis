#include "SpectrumAnalyzer.h"

// print all adjacent vertices of given vertex
template <typename T>
void display_AdjList(adjNode<T>* ptr, int i)
{
    while (ptr != nullptr) {
        std::cout << "(" << i << ", " << ptr->val
                  << ", " << ptr->cost << " , " << ptr->payload <<") ";
        ptr = ptr->next;
    }
    std::cout << std::endl;
}

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name) {
    TFile * file = new TFile(file_name.c_str());
    TH2D gg = *((TH2D*) file->Get("mgamma_gamma"));

    file->Close();

    double levels [] = {0, 121.8, 366.5, 1085.5, 1233.9, 1529.8};
    int gammas [][6] = { {5, 1, 3},
                        {4, 1},
                        {3, 0, 1},
                        {2, 1},
                        {1, 0},
                        {0}};
    graphEdge<double> edges[] = {
            // (x, y, w) -> edge from x to y with weight w
            {5,1,0},
            {5,3,0},
            {4,1,0},
            {3,0,0},
            {3,1,0},
            {2,1,0},
            {1,0,0},
            {0,0,0}
    };

    DiaGraph<double> test(edges, levels, [](double a, double b){return a-b;}, 8,6);
    for (unsigned ii=0; ii<6; ++ii)
        display_AdjList(test.head[ii], ii);
}