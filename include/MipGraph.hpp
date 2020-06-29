#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <ostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <NodeBddTable.hpp>
#include <NodeBdd.hpp>
#include <vector>
#include "ZddNode.hpp"
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

typedef property < vertex_index_t, int,
        property < vertex_name_t, NodeId,
        property < vertex_degree_t, int,
        property < vertex_distance_t, VarsNode,
        property < vertex_color_t, std::shared_ptr<SubNodeZdd<>> 
        > > > > > VertexProperty;

typedef property < edge_index_t, int, 
        property < edge_weight_t, bool,
        property < edge_weight2_t, VarsEdge
        > > > EdgeProperty;

using MipGraph = adjacency_list<vecS,vecS,bidirectionalS,VertexProperty,EdgeProperty>;
using Edge = graph_traits<MipGraph>::edge_descriptor;
using Vertex = graph_traits<MipGraph>::vertex_descriptor;

typedef property_map<MipGraph, vertex_index_t>::type IndexAccessor;
typedef property_map<MipGraph, vertex_name_t>::type NodeIdAccessor;
typedef property_map<MipGraph, vertex_color_t>::type NodeZddIdAccessor;
typedef property_map<MipGraph, vertex_degree_t>::type NodeMipIdAccessor;
typedef property_map<MipGraph, vertex_distance_t>::type VarsNodeAccessor;
typedef property_map<MipGraph, edge_index_t>::type EdgeIndexAccessor;
typedef property_map<MipGraph, edge_weight_t>::type EdgeTypeAccessor;
typedef property_map<MipGraph, edge_weight2_t>::type EdgeVarAccessor;
typedef std::vector<std::vector<double>> dbl_matrix;
class ColorWriterEdgeX {
private:
    const MipGraph& g;
    const NodeTableEntity<>* table;
    const dbl_matrix* lp_x_high;
    const dbl_matrix* lp_x_low;

public:
    explicit ColorWriterEdgeX(MipGraph& _g, NodeTableEntity<>* _table, dbl_matrix* _lp_high, dbl_matrix* _lp_low) : g{_g},table(_table), lp_x_high(_lp_high), lp_x_low(_lp_low) {

    }

    void operator()(std::ostream &output, Edge _edge) {
        auto high =  get(boost::edge_weight_t(),g,_edge);
        auto node_id = get(boost::vertex_name_t(),g,source(_edge,g));
        auto node = table->node(node_id);
        int t = node.get_weight();
        int j = node.get_job()->job;
        if(high) {
            if ((*lp_x_high)[j][t] > 0.00001) {
                output << "[label = "<< (*lp_x_high)[j][t] << ",color = red]";
            } else {
                output << "[label = "<< (*lp_x_high)[j][t] <<"]";
            }
        } else  {
            if((*lp_x_low)[j][t] > 0.00001) {
                output << "[label = "<< (*lp_x_low)[j][t] <<",color = red, style = dashed]";
            } else {
                output << "[label = "<< (*lp_x_low)[j][t] <<",style=dashed]";
            }
        }     
    }
};

class ColorWriterEdgeIndex {
private:
    const MipGraph& g;

public:
    explicit ColorWriterEdgeIndex(MipGraph& _g) : g{_g} {

    }

    void operator()(std::ostream &output, Edge _edge) {
        int index = get(boost::edge_index_t(),g,_edge);
        bool high =  get(boost::edge_weight_t(),g,_edge);       
    
            if(high) {
                output << "[label = "<< index << "]";
            } else {
                output << "[label = "<< index <<", style = dashed]";
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
        NodeId tmp_nodeid = get(boost::vertex_name_t(), g, _vertex);
        if(tmp_nodeid > 1) {
            output <<  " [label=\" " << table.node(tmp_nodeid).get_job()->job << " " << table.node(tmp_nodeid).get_weight()  << "\"]";
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
