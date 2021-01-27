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

// using namespace boost;

struct VarsEdge {
    GRBVar x;
    GRBVar alpha;
};

struct VarsNode {
    std::array<GRBVar, 2> omega;
};

using VertexProperty = boost::property<
    boost::vertex_index_t,
    int,
    boost::property<
        boost::vertex_name_t,
        NodeId,
        boost::property<
            boost::vertex_degree_t,

            int,
            boost::property<boost::vertex_distance_t,
                            VarsNode,
                            boost::property<boost::vertex_color_t,
                                            std::shared_ptr<SubNodeZdd<>>>>>>>;

using EdgeProperty = boost::property<
    boost::edge_index_t,
    int,
    boost::property<boost::edge_weight_t,
                    bool,
                    boost::property<boost::edge_weight2_t, VarsEdge>>>;

using MipGraph = boost::adjacency_list<boost::vecS,
                                       boost::vecS,
                                       boost::bidirectionalS,
                                       VertexProperty,
                                       EdgeProperty>;
using Edge = boost::graph_traits<MipGraph>::edge_descriptor;
using Vertex = boost::graph_traits<MipGraph>::vertex_descriptor;

using IndexAccessor =
    boost::property_map<MipGraph, boost::vertex_index_t>::type;
using NodeIdAccessor =
    boost::property_map<MipGraph, boost::vertex_name_t>::type;
using NodeZddIdAccessor =
    boost::property_map<MipGraph, boost::vertex_color_t>::type;
using NodeMipIdAccessor =
    boost::property_map<MipGraph, boost::vertex_degree_t>::type;
using VarsNodeAccessor =
    boost::property_map<MipGraph, boost::vertex_distance_t>::type;
using EdgeIndexAccessor =
    boost::property_map<MipGraph, boost::edge_index_t>::type;
using EdgeTypeAccessor =
    boost::property_map<MipGraph, boost::edge_weight_t>::type;
using EdgeVarAccessor =
    boost::property_map<MipGraph, boost::edge_weight2_t>::type;
using dbl_matrix = std::vector<std::vector<double>>;
class ColorWriterEdgeX {
   private:
    const MipGraph&          g;
    const NodeTableEntity<>* table;
    static constexpr double  EPS_GRAPH = 1e-6;

   public:
    explicit ColorWriterEdgeX(const MipGraph&          _g,
                              const NodeTableEntity<>* _table)
        : g{_g},
          table(_table) {}

    void operator()(std::ostream& output, Edge _edge) {
        auto  high = get(boost::edge_weight_t(), g, _edge);
        auto  node_id = get(boost::vertex_name_t(), g, source(_edge, g));
        auto& node = table->node(node_id);

        if (high) {
            auto& x = node.lp_x[1];
            if (x > EPS_GRAPH) {
                output << "[label = " << x << ",color = red]";
            } else {
                output << "[label = " << x << "]";
            }
        } else {
            auto& x = node.lp_x[0];
            if (x > EPS_GRAPH) {
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
    explicit ColorWriterEdgeIndex(const MipGraph& _g) : g{_g} {}

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
    const MipGraph&          g;
    const NodeTableEntity<>& table;

   public:
    ColorWriterVertex(const MipGraph& _g, const NodeTableEntity<>& _table)
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
