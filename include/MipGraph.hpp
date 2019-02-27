#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <ostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <NodeBddTable.hpp>
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

typedef property_map<MipGraph, vertex_index_t>::type IndexAccessor;
typedef property_map<MipGraph, vertex_name_t>::type NodeIdAccessor;
typedef property_map<MipGraph, edge_index_t>::type EdgeIndexAccessor;
typedef property_map<MipGraph, edge_weight_t>::type EdgeTypeAccessor;
typedef property_map<MipGraph, edge_weight2_t>::type EdgeVarAccessor;

class ColorWriterEdge {
private:
    const MipGraph& g;
    const double *x;

public:
    explicit ColorWriterEdge(MipGraph& _g, double *_x) : g{_g}, x{_x} {

    }

    void operator()(std::ostream &output, Edge _edge) {       
        if(x[get(boost::edge_index_t(),g,_edge)] > 0.00001) {
            if(get(boost::edge_index_t(),g,_edge)%2 == 0) {
                output << "[label = "<< x[get(boost::edge_index_t(),g,_edge)] <<",color = red]";
            } else {
                output << "[label = "<< x[get(boost::edge_index_t(),g,_edge)] <<",color = red, style = dashed]";
            }
        } else {
            output << "[label = "<< x[get(boost::edge_index_t(),g,_edge)] <<"]";
        }
    }
};

class ColorWriterVertex {
private:
    const MipGraph &g;
    NodeTableEntity<>& table;
public:
    ColorWriterVertex(MipGraph &_g, NodeTableEntity<>& _table) : g{_g}, table{_table} {

    }

    void operator()(std::ostream &output, Vertex _vertex) {
        nodeid tmp_nodeid = get(boost::vertex_name_t(), g, _vertex);
        if(tmp_nodeid > 1) {
            output <<  " [label=\" " << table.node(tmp_nodeid).GetJob()->job << " " << table.node(tmp_nodeid).GetWeight()  << "\"]";
        }
    }
};

#endif // MIP_GRAPH_HPP
