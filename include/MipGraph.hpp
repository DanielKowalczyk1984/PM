#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <gurobi_c++.h>
#include <NodeBddTable.hpp>
#include <NodeId.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <ostream>

struct VertexData {
    size_t index{};
    NodeId node_id{};
};

struct EdgeData {
    size_t id{};
    bool   high{};
    GRBVar x{};
};
using MipGraph = boost::adjacency_list<boost::vecS,
                                       boost::vecS,
                                       boost::bidirectionalS,
                                       VertexData,
                                       EdgeData>;
using Edge = boost::graph_traits<MipGraph>::edge_descriptor;
using Vertex = boost::graph_traits<MipGraph>::vertex_descriptor;

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

    void operator()(std::ostream& output, const Edge& _edge) {
        auto  node_id = g[source(_edge, g)].node_id;
        auto& node = table->node(node_id);

        if (g[_edge].high) {
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

    void operator()(std::ostream& output, const Edge& _edge) {
        auto index = g[_edge].id;
        auto high = g[_edge].high;

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

    void operator()(std::ostream& output, const Vertex& _vertex) {
        NodeId tmp_nodeid = g[_vertex].node_id;
        if (tmp_nodeid > 1) {
            output << " [label=\" " << table.node(tmp_nodeid).get_job()->job
                   << " " << table.node(tmp_nodeid).get_weight() << "\"]";
        }
    }
};

#endif  // MIP_GRAPH_HPP
