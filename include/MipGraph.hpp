#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <ostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <NodeBddTable.hpp>
#include <node_duration.hpp>
#include <gurobi_c++.h>
#include <scheduleset.h>

using namespace boost;

struct VarsEdge {
    GRBVar x;
    GRBVar alpha;
};

struct VarsNode{
    GRBVar omega[2];
};

typedef adjacency_list_traits < vecS, vecS, bidirectionalS > Traits;

typedef property < vertex_index_t, int,
        property < vertex_name_t, nodeid,
        property < vertex_distance_t, VarsNode
        > > > VertexProperty;

typedef property < edge_index_t, int, 
        property < edge_weight_t, bool,
        property < edge_weight2_t, VarsEdge
        > > > EdgeProperty;

using MipGraph = adjacency_list<vecS,vecS,bidirectionalS,VertexProperty,EdgeProperty>;
using Edge = graph_traits<MipGraph>::edge_descriptor;
using Vertex = graph_traits<MipGraph>::vertex_descriptor;
using EdgeIterator = graph_traits<MipGraph>::edge_iterator;

typedef property_map<MipGraph, vertex_index_t>::type IndexAccessor;
typedef property_map<MipGraph, vertex_name_t>::type NodeIdAccessor;
typedef property_map<MipGraph, vertex_distance_t>::type VarsNodeAccessor;
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
        int index = get(boost::edge_index_t(),g,_edge);
        bool high =  get(boost::edge_weight_t(),g,_edge);       
        if(x[index] > 0.00001) {
            if(high) {
                output << "[label = "<< index << ",color = red, style =dashed]";
            } else {
                output << "[label = "<< index <<",color = red]";
            }
        } else {
            if(high) {
                output << "[label = "<< index <<",style = dashed]";
            } else {
                output << "[label = "<< index <<"]";
            }
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

/*class ColorWriterSolution {
private:
    const MipGraph& g;

public:
    explicit ColorWriterSolution(MipGraph& _g) : g{_g} {

    }

    void operator()(std::ostream &output, Edge _edge) {
        VarsEdge var =  get(boost::edge_weight2_t(),g,_edge);
        bool high =  get(boost::edge_weight_t(),g,_edge);
        double x = var.x.get(GRB_DoubleAttr_X);  
        if(x > 0.00001) {
            if(high) {
                output << "[label = "<< x << ",color = red, style =dashed]";
            } else {
                output << "[label = "<< x <<",color = red]";
            }
        } else {
            if(high) {
                output << "[label = "<< x <<",style = dashed]";
            } else {
                output << "[label = "<< x <<"]";
            }
        }
    }
};*/

#endif // MIP_GRAPH_HPP
