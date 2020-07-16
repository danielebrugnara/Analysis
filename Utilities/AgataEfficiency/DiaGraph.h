#ifndef EFFANALYSIS_DIAGRAPH_H
#define EFFANALYSIS_DIAGRAPH_H

#include <iostream>
#include <vector>
// stores adjacency list items
template <typename Tw>
struct adjEdge {
    int pointed_ver;
    Tw weight;
    adjEdge<Tw>* next;
    friend std::ostream& operator<<(std::ostream& os, const adjEdge<Tw>& ed){
        os << "(edg pts. to : " << ed.pointed_ver << " , contains : {"<< ed.weight <<"} )";
        return os;
    }

};

template <typename Tw, typename Tp>
struct graphNode{
    Tp payload;
    adjEdge<Tw>* adj;
    friend std::ostream& operator<<(std::ostream& os, const graphNode<Tw, Tp>& nd){
        os << "---->(node contains : {"<< nd.payload <<"} )-->";
        return os;
    }
};

// structure to store edges
template <typename Tw>
struct graphEdge {
    int start_ver, end_ver;
    Tw weight;
};

template <typename Tw, typename  Tp>
class DiaGraph{
public:
    DiaGraph(graphEdge<Tw> edges[], Tp* payloads, int n, int N):N(N){
        // allocate new node
        head = new graphNode<Tw, Tp>[N];
        // initialize head pointer for all vertices
        for (int i = 0; i < N; ++i)
            head[i] = {payloads[i],nullptr};
        // construct directed graph by adding edges to it
        for (int i = 0; i < n; ++i)  {
            int start_ver = edges[i].start_ver;
            int end_ver = edges[i].end_ver;
            Tw weight = edges[i].weight;
            // insert in the beginning
            head[start_ver].adj = addAdjListNode(end_ver, weight, head[start_ver].adj);
        }
    }
    // Destructor
    ~DiaGraph() {
        for (int i = 0; i < N; i++){
            adjEdge<Tw>* ptr = head[i].adj;
            adjEdge<Tw>* ptr_next;
            while (ptr != nullptr) {
                ptr_next = ptr->next;
                delete ptr;
                ptr = ptr_next;
            }
        }
        delete[] head;
    }

    int GetNumberNodes(){return N;}

    std::vector<std::pair<const adjEdge<Tw>*,const adjEdge<Tw>*>> getConsecutiveEdges(const int& startNodeidx){
        std::vector<std::pair<const adjEdge<Tw>*,const adjEdge<Tw>*>> list;
        auto* gamma1 = head[startNodeidx].adj;
        while(gamma1 != nullptr){
            auto* gamma2 = head[gamma1->pointed_ver].adj;
            while(gamma2 != nullptr){
                list.emplace_back(gamma1, gamma2);
                gamma2 = gamma2->next;
            }
            gamma1 = gamma1->next;
        }
        return list;
    }
private:
    graphNode<Tw, Tp> *head;
    int N;  // number of nodes in the graph
    // insert new nodes into adjacency list from given graph
    adjEdge<Tw>* addAdjListNode(int pointed_ver, Tw weight, adjEdge<Tw>* head)   {
        return new adjEdge<Tw>({pointed_ver, weight, head});
    }


    friend std::ostream& operator<<(std::ostream& os, const DiaGraph<Tw, Tp>& gr)
    {
        for(int ii=0; ii<gr.N; ++ii) {
            adjEdge<Tw>* ptr = gr.head[ii].adj;
            os <<"idx. " << ii <<" : "<<gr.head[ii];
            while (ptr != nullptr) {
                os<< *ptr << "-->";
                ptr = ptr->next;
            }
            os << " nullptr "<< std::endl;
        }
        return os;
    }
};

#endif //EFFANALYSIS_DIAGRAPH_H
