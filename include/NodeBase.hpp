#ifndef NODE_BASE_HPP
#define NODE_BASE_HPP

#include "OptimalSolution.hpp"
#include "NodeBdd.hpp"

using NodeId = tdzdd::NodeId;

// template<typename T = double>
class NodeBase {
  private:
    int num_layer;

    bool root_node;
    bool terminal_node;

    Job *job;

  public:
    NodeId branch[2];

    /**
     * Constructor
     */
    NodeBase():
          num_layer(0),
          root_node(false),
          terminal_node(false),
          job(nullptr),
          branch{NodeId(),NodeId()}
          {

          };

    NodeBase(int &_num_layer, bool &_root_node,bool &_terminal_node):
         num_layer(_num_layer),
         root_node(_root_node),
         terminal_node(_terminal_node),
         job(nullptr),
         branch{NodeId(),NodeId()} {
    }

    NodeBase(int i, int j) :
          num_layer(0),
          root_node(false),
          terminal_node(false),
          job(nullptr),
          branch{i,j}{
      }

    NodeBase(const NodeBase &src) = default;
    NodeBase(NodeBase &&src) = default;
    NodeBase& operator=(const NodeBase &src) = default;
    NodeBase& operator=(NodeBase &&src) = default;

    void set_job(Job *_job, bool _terminal_node = false){
      job = _job;
      terminal_node = _terminal_node;
    }

    void set_layer(int _num_layer) {
      num_layer = _num_layer;
    }

    void set_root_node (bool _root_node) {
      root_node = _root_node;
    }

    void set_terminal_node (bool _terminal_node) {
      root_node = _terminal_node;
    }

    int get_layer_number() {
      return num_layer;
    }

    bool is_terminal_node(){
      return terminal_node;
    }

    Job *get_job(){
      return job;
    }

    size_t hash() const {
        size_t h = branch[0].code();
        for (int i = 1; i < 2; ++i) {
            h = h * 314159257 + branch[i].code() * 271828171;
        }
        return h;
    }

    bool operator==(NodeBase const& o) const {
        for (int i = 0; i < 2; ++i) {
            if (branch[i] != o.branch[i]) return false;
        }
        return true;
    }

    bool operator!=(NodeBase const& o) const {
        return !operator==(o);
    }

    friend std::ostream& operator<<(std::ostream& os, NodeBase const& o) {
        os << "(" << o.branch[0];
        for (int i = 1; i < 2; ++i) {
            os << "," << o.branch[i];
        }
        return os << ")";
    }

};


#endif // NODE_BASE_HPP



