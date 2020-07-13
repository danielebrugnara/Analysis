#ifndef EFFANALYSIS_DIAGRAPH_H
#define EFFANALYSIS_DIAGRAPH_H

#include <iostream>
// stores adjacency list items
template <typename T>
struct adjNode {
    int val;
    T cost;
    T payload;
    adjNode* next;
};

// structure to store edges
template <typename T>
struct graphEdge {
    int start_ver, end_ver;
    T weight;
};

template <typename T>
class DiaGraph{
    // insert new nodes into adjacency list from given graph
    adjNode<T>* getAdjListNode(int value, int weight, adjNode<T>* head)   {
        adjNode<T>* newNode = new adjNode<T>;
        newNode->val = value;
        newNode->cost = weight;
        newNode->payload = weight;

        newNode->next = head;   // point new node to current head
        return newNode;
    }
    adjNode<T>* getAdjListNode(int value, T weight, T payload, adjNode<T>* head)   {
        adjNode<T>* newNode = new adjNode<T>;
        newNode->val = value;
        newNode->cost = weight;
        newNode->payload = payload;

        newNode->next = head;   // point new node to current head
        return newNode;
    }
    int N;  // number of nodes in the graph
public:
    adjNode<T> **head;                //adjacency list as array of pointers
    // Constructor
    DiaGraph(graphEdge<T> edges[], int n, int N)  {
        // allocate new node
        head = new adjNode<T>*[N]();
        this->N = N;
        // initialize head pointer for all vertices
        for (int i = 0; i < N; ++i)
            head[i] = nullptr;
        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < n; ++i)  {
            int start_ver = edges[i].start_ver;
            int end_ver = edges[i].end_ver;
            int weight = edges[i].weight;
            // insert in the beginning
            adjNode<T>* newNode = getAdjListNode(end_ver, weight, head[start_ver]);

            // point head pointer to new node
            head[start_ver] = newNode;
        }
    }

    DiaGraph(graphEdge<T> edges[], T* payloads, T compute_weight(T, T), int n, int N){
        // allocate new node
        head = new adjNode<T>*[N]();
        this->N = N;
        // initialize head pointer for all vertices
        for (int i = 0; i < N; ++i)
            head[i] = nullptr;
        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < n; ++i)  {
            int start_ver = edges[i].start_ver;
            int end_ver = edges[i].end_ver;
            T weight = edges[i].weight;
            // insert in the beginning
            adjNode<T>* newNode = getAdjListNode(end_ver, weight, payloads[start_ver], head[start_ver]);

            // point head pointer to new node
            head[start_ver] = newNode;
        }
        for (unsigned i = 0; i < N; ++i)  {
            adjNode<T>* ptr=head[i];
            while (ptr != nullptr) {
                ptr->cost = compute_weight(head[i]->payload, ptr->payload);
                ptr = ptr->next;
            }

        }

    }
    // Destructor
    ~DiaGraph() {
        for (int i = 0; i < N; i++)
            delete[] head[i];
        delete[] head;
    }
};

#endif //EFFANALYSIS_DIAGRAPH_H
