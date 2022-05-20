// graph.h <Starter Code>
// Harith Patel
//
// Basic graph class using adjacency list representation.  
//
// University of Illinois at Chicago
// CS 251: Spring 2022
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <map>
#include <vector>
#include <set>
#include <utility>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
  map<VertexT, map<VertexT, WeightT> > AdjList;

 public:
  //
  // defult constructor:
  //
  // Constructs an empty graph where
  //
  graph() = default;

  //
  // NumVertices:
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return static_cast<int>(AdjList.size());
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    int count = 0;

    // loop through the adjacency list and count how many
    // edges currently exist
    for (auto ver : AdjList)
      count += ver.second.size();

    return count;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    if (AdjList.count(v) == 1) return false;

    // empty map for edges
    map<VertexT, WeightT> edges;

    // if we get here, vertex does not exist so insert.

    AdjList.insert(pair<VertexT, map<VertexT, WeightT> >(v, edges));

    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    auto itr = AdjList.find(from);
    if (itr == AdjList.end()) return false;
    if (AdjList.count(to) != 1) return false;

    if (itr->second.count(to) == 1)
      itr->second[to] = weight;
    else
      itr->second.insert(pair<VertexT, WeightT>(to, weight));

    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    auto itr = AdjList.find(from);
    if (itr == AdjList.end()) return false;
    auto itr2 = itr->second.find(to);
    if (itr2 == itr->second.end()) return false;
    weight = itr2->second;

    return true;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT>  S;

    auto itr = AdjList.find(v);
    if (itr == AdjList.end()) return S;

    for (auto v : itr->second) 
      S.insert(v.first);

    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> v;

    for (auto ver : AdjList)
      v.push_back(ver.first);

    return v;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:Edges:**" << endl;
    vector<VertexT> vert = getVertices();

    for (auto v : vert) {
      cout << v << ": ";
      set<VertexT> edges = neighbors(v);

      for (auto e : edges) {
        WeightT weight;

        getWeight(v, e, weight);
        cout << "(" << v << "," << e << "," << weight << ") ";
      }

      cout << endl;
    }

    output << "**************************************************" << endl;
  }
};
