#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <gurobi_c++.h>
#include <scheduleset.h>
#include <NodeBdd.hpp>
#include <NodeBddTable.hpp>
#include <array>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_set.hpp>
#include <ostream>
#include <vector>
#include "ZddNode.hpp"

using namespace boost;

struct VarsEdge {
    GRBVar x;
    GRBVar alpha;
};

struct VarsNode {
    std::array<GRBVar, 2> omega;
};

using VertexProperty = property<
    vertex_index_t,
    int,
    property<vertex_name_t,
             NodeId,
             property<vertex_degree_t,
                      int,
                      property<vertex_distance_t,
                               VarsNode,
                               property<vertex_color_t,
                                        std::shared_ptr<SubNodeZdd<>>>>>>>;

using EdgeProperty =
    property<edge_index_t,
             int,
             property<edge_weight_t, bool, property<edge_weight2_t, VarsEdge>>>;

using MipGraph =
    adjacency_list<vecS, vecS, bidirectionalS, VertexProperty, EdgeProperty>;
using Edge = graph_traits<MipGraph>::edge_descriptor;
using Vertex = graph_traits<MipGraph>::vertex_descriptor;

using IndexAccessor = property_map<MipGraph, vertex_index_t>::type;
using NodeIdAccessor = property_map<MipGraph, vertex_name_t>::type;
using NodeZddIdAccessor = property_map<MipGraph, vertex_color_t>::type;
using NodeMipIdAccessor = property_map<MipGraph, vertex_degree_t>::type;
using VarsNodeAccessor = property_map<MipGraph, vertex_distance_t>::type;
using EdgeIndexAccessor = property_map<MipGraph, edge_index_t>::type;
using EdgeTypeAccessor = property_map<MipGraph, edge_weight_t>::type;
using EdgeVarAccessor = property_map<MipGraph, edge_weight2_t>::type;
using dbl_matrix = std::vector<std::vector<double>>;
class ColorWriterEdgeX {
   private:
    const MipGraph&          g;
    const NodeTableEntity<>* table;

   public:
    explicit ColorWriterEdgeX(MipGraph& _g, NodeTableEntity<>* _table)
        : g{_g},
          table(_table) {}

    void operator()(std::ostream& output, Edge _edge) {
        auto  high = get(boost::edge_weight_t(), g, _edge);
        auto  node_id = get(boost::vertex_name_t(), g, source(_edge, g));
        auto& node = table->node(node_id);
        auto& x = node.lp_x[high];

        if (high) {
            if (x > 1e-5) {
                output << "[label = " << x << ",color = red]";
            } else {
                output << "[label = " << x << "]";
            }
        } else {
            if (x > 1e-5) {
                output << "[label = " << x << ",color = red, style = dashed]";
            } else {
                output << "[label = " << x << ",style=dashed]";
            }
        }
    }
};

class ColorWriterEdgeIndex {
   private:
    const MipGraph& g;

   public:
    explicit ColorWriterEdgeIndex(MipGraph& _g) : g{_g} {}

    void operator()(std::ostream& output, Edge _edge) {
        int  index = get(boost::edge_index_t(), g, _edge);
        bool high = get(boost::edge_weight_t(), g, _edge);

        if (high) {
            output << "[label = " << index << "]";
        } else {
            output << "[label = " << index << ", style = dashed]";
        }
    }
};

class ColorWriterVertex {
   private:
    const MipGraph&    g;
    NodeTableEntity<>& table;

   public:
    ColorWriterVertex(MipGraph& _g, NodeTableEntity<>& _table)
        : g{_g},
          table{_table} {}

    void operator()(std::ostream& output, Vertex _vertex) {
        NodeId tmp_nodeid = get(boost::vertex_name_t(), g, _vertex);
        if (tmp_nodeid > 1) {
            output << " [label=\" " << table.node(tmp_nodeid).get_job()->job
                   << " " << table.node(tmp_nodeid).get_weight() << "\"]";
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

#endif  // MIP_GRAPH_HPP
