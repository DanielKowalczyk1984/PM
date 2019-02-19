#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <node_duration.hpp>
#include <gurobi_c++.h>

using namespace boost;

typedef adjacency_list_traits < vecS, vecS, bidirectionalS > Traits;

typedef property < vertex_index_t, int,
        property < vertex_name_t, nodeid
        > > VertexProperty;

typedef property < edge_index_t, int, 
        property < edge_weight_t, bool,
        property < edge_weight2_t, GRBVar
        > > > EdgeProperty;

using MipGraph = adjacency_list<vecS,vecS,bidirectionalS,VertexProperty,EdgeProperty>;
using Edge = graph_traits<MipGraph>::edge_descriptor;
using Vertex = graph_traits<MipGraph>::vertex_descriptor;
using EdgeIterator = graph_traits<MipGraph>::edge_iterator;

typedef property_map<MipGraph, vertex_name_t>::type NodeIdAccessor;
typedef property_map<MipGraph, edge_weight_t>::type EdgeTypeAccessor;
typedef property_map<MipGraph, edge_weight2_t>::type EdgeVarAccessor;

#endif // MIP_GRAPH_HPP
