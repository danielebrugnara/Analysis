#ifndef EFFANALYSIS_DIAGRAPH_H
#define EFFANALYSIS_DIAGRAPH_H

#include <iostream>
// stores adjacency list items
template <typename Tw>
struct adjEdge {
    int pointed_ver;
    Tw weight;
    adjEdge<Tw>* next;
};

template <typename Tw, typename Tp>
struct graphNode{
    Tp payload;
    adjEdge<Tw>* adj;
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
    //adjEdge<Tw> **head;                //adjacency list as array of pointers
    graphNode<Tw, Tp> *head;
private:
    int N;  // number of nodes in the graph
    // insert new nodes into adjacency list from given graph
    adjEdge<Tw>* addAdjListNode(int pointed_ver, Tw weight, adjEdge<Tw>* head)   {
        return new adjEdge<Tw>({pointed_ver, weight, head});
    }

public:
    friend std::ostream& operator<<(std::ostream& os, const DiaGraph<Tw, Tp>& gr)
    {
        for(int ii=0; ii<gr.N; ++ii) {
            adjEdge<Tw>* ptr = gr.head[ii].adj;
            while (ptr != nullptr) {
                os << "(" << ii << ", " << ptr->pointed_ver
                          << ", " << ptr->weight << " , " << gr.head[ii].payload << ") ";
                ptr = ptr->next;
            }
            os << std::endl;
        }
        return os;
    }
};

#endif //EFFANALYSIS_DIAGRAPH_H
