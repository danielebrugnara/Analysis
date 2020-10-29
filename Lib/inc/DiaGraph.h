#pragma once

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

//Tw stants for weight (link) content, Tp stands for (node) payload content
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
        DiaGraph(std::vector<graphEdge<Tw>> edges, std::vector<Tp> payloads):N(payloads.size()){
            // allocate new node
            head = new graphNode<Tw, Tp>[N];
            // initialize head pointer for all vertices
            for (long unsigned int i = 0; i < payloads.size(); ++i)
                head[i] = {payloads[i],nullptr};
            // construct directed graph by adding edges to it
            for (long unsigned int i = 0; i < edges.size(); ++i)  {
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

        std::vector<std::pair<const Tw*,const Tw*>> getConsecutiveEdgeWeights(const int& startNodeidx){
            std::vector<std::pair<const Tw*,const Tw*>> list;
            auto* edge1 = head[startNodeidx].adj;
            while(edge1 != nullptr){
                auto* edge2 = head[edge1->pointed_ver].adj;
                while(edge2 != nullptr){
                    list.emplace_back(&edge1->weight, &edge2->weight);
                    edge2 = edge2->next;
                }
                edge1 = edge1->next;
            }
            return list;
        }

        std::vector<std::pair<const Tw*, const Tp*>> getAdjacent(const int& Nodeidx){
            std::vector<std::pair<const Tw*, const Tp*>> list;
            auto* edge = head[Nodeidx].adj;
            while(edge != nullptr){
                list.emplace_back(&edge->weight, &head[edge->pointed_ver].payload);
                edge = edge->next;
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
    public:
////        template<typename Tw>
//            class weight_iterator
//            {
//                public:
//
//                    using iterator_category = std::random_access_iterator_tag;
//                    using value_type = adjEdge<Tw>;
//                    using difference_type = std::ptrdiff_t;
//                    using pointer = adjEdge<Tw>*;
//                    using reference = adjEdge<Tw>&;
//
//                public:
//
//                    weight_iterator(Tw* ptr = nullptr){m_ptr = ptr;}
//                    weight_iterator(const weight_iterator<Tw>& rawIterator) = default;
//                    ~weight_iterator(){}
//
//                    weight_iterator<Tw>&                  operator=(const weight_iterator<Tw>& rawIterator) = default;
//                    weight_iterator<Tw>&                  operator=(Tw* ptr){m_ptr = ptr;return (*this);}
//
//                    operator                                    bool()const
//                    {
//                        if(m_ptr)
//                            return true;
//                        else
//                            return false;
//                    }
//
//                    bool                                        operator==(const weight_iterator<Tw>& rawIterator)const{return (m_ptr == rawIterator.getConstPtr());}
//                    bool                                        operator!=(const weight_iterator<Tw>& rawIterator)const{return (m_ptr != rawIterator.getConstPtr());}
//
//                    weight_iterator<Tw>&                  operator+=(const difference_type& movement){for(int i=0; i<movement;++i)
//                                                                                                        m_ptr += movement;return (*this);
//                                                                                                     }
//                    weight_iterator<Tw>&                  operator-=(const difference_type& movement){m_ptr -= movement;return (*this);}
//                    weight_iterator<Tw>&                  operator++(){++m_ptr;return (*this);}
//                    weight_iterator<Tw>&                  operator--(){--m_ptr;return (*this);}
//                    weight_iterator<Tw>                   operator++(int){auto temp(*this);++m_ptr;return temp;}
//                    weight_iterator<Tw>                   operator--(int){auto temp(*this);--m_ptr;return temp;}
//                    weight_iterator<Tw>                   operator+(const difference_type& movement){auto oldPtr = m_ptr;m_ptr+=movement;auto temp(*this);m_ptr = oldPtr;return temp;}
//                    weight_iterator<Tw>                   operator-(const difference_type& movement){auto oldPtr = m_ptr;m_ptr-=movement;auto temp(*this);m_ptr = oldPtr;return temp;}
//
//                    difference_type                     operator-(const weight_iterator<Tw>& rawIterator){return std::distance(rawIterator.getPtr(),this->getPtr());}
//
//                    Tw&                                 operator*(){return *m_ptr;}
//                    const Tw&                           operator*()const{return *m_ptr;}
//                    Tw*                                 operator->(){return m_ptr;}
//
//                    Tw*                                 getPtr()const{return m_ptr;}
//                    const Tw*                           getConstPtr()const{return m_ptr;}
//
//                protected:
//
//                    Tw*                                 m_ptr;
//            };
};