#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <memory>
#include <OptimalSolution.hpp>

template<typename T>
class NodeDuration;

template<typename T>
class PrevNode {
  private:
    PrevNode<T> *prev;
    bool high;
    NodeDuration<T>* head_node;

  public:
    T f;
    /**
     * Constructor
     */
    PrevNode(T &_f, PrevNode<T> *&_prev, bool &_high): f(_f),
        prev(_prev),
        high(_high){

    };
    PrevNode(){

    };

    /**
     * Copy Constructor
     */
    PrevNode<T>(PrevNode<T> &src): f(src.f), prev(src.prev), high(src.high), head_node(src.head_node){
      
    }

    /**
     * Move Constructor
     */
    PrevNode<T>(PrevNode<T> &&src): f(src.f), prev(src.prev), high (src.high), head_node (src.head_node){
      
    }

    /**
     * Copy Assignment
     */
    PrevNode<T>& operator=(const PrevNode<T> &src){
      if(&src == this) {
        return *this;
      }

      f = src.f;
      prev = src.prev;
      high = src.high;
      head_node = src.head_node;

      return *this;
    }

    /**
     * Move Assignment
     */
    PrevNode<T>& operator=(const PrevNode<T> &&src){
      if(&src == this) {
        return *this;
      }

      f = src.f;
      prev = src.prev;
      high = src.high;
      head_node = src.node;

      return *this;
    }

    void SetPrev(PrevNode<T> * _prev) {
        prev = _prev;
    }

    void SetF(T _f) {
        f = _f;
    }

    void SetHigh(bool _high) {
        high = _high;
    }

    void SetHeadNode(NodeDuration<T>* _head){
      head_node = _head;
    }

    void Reset(){
      f = -DBL_MAX;
      prev = nullptr;
      high = false;
    }

    T GetF() const {
        return f;
    }

    PrevNode<T>* GetPrev() {
        return prev;
    }

    bool& GetHigh() {
        return high;
    }

    Job* GetJob() {
        return head_node->GetJob();
    }

    int& GetWeight(){
      return head_node->GetWeight();
    }
};

template<typename T>
class NodeDuration {
  private:
    int weight;
    int num_layer;

    bool root_node;
    bool terminal_node;

    Job *job;

  public:
    PrevNode<T> prev1;
    PrevNode<T> prev2;

    NodeDuration<T>* y;
    NodeDuration<T>* n;
    
    /**
     * Constructor
     */
    NodeDuration()
        : num_layer(0),
          weight(0),
          root_node(false),
          terminal_node(false) {

    };

    NodeDuration(int &_weight, int &_num_layer, bool &_root_node,bool &_terminal_node):
         weight(_weight), num_layer(_num_layer), root_node(_root_node),terminal_node(_terminal_node) {
          prev1.Reset();
          prev2.Reset();
    }

    /**
     * Copy Constructor
     */
    NodeDuration<T>(const NodeDuration<T> &src) {
      weight = src.weight;
      num_layer = src.num_layer;
      root_node = src.root_node;
      terminal_node = src.terminal_node;
      job = src.job;

      y = src.y;
      n = src.n;

      prev1 = src.prev1;
      prev2 = src.prev2;
    }

    /**
     * Move Constructor
     */
    NodeDuration<T>(NodeDuration<T> &&src) {
      weight = src.weight;
      num_layer = src.num_layer;
      root_node = src.root_node;
      terminal_node = src.terminal_node;
      job = src.job;

      y = src.y;
      n = src.n;

      prev1 = src.prev1;
      prev2 = src.prev2;
    }

    /**
     * Copy Operator
     */
    NodeDuration<T>& operator=(const NodeDuration<T> &src){
      if(&src == this) {
        return *this;
      }

      weight = src.weight;
      num_layer = src.num_layer;
      root_node = src.root_node;
      terminal_node = src.node;
      job = src.job;

      y = src.y;
      n = src.n;

      prev1 = src.prev1;
      prev2 = src.prev2;

      return *this;
    }

    /**
     * Move Operator
     */
    NodeDuration<T>& operator=(NodeDuration<T> &&src){
      if(&src == this) {
        return *this;
      }

      weight = src.weight;
      num_layer = src.num_layer;
      root_node = src.root_node;
      terminal_node = src.node;
      job = src.job;

      y = src.y;
      n = src.n;

      prev1 = src.prev1;
      prev2 = src.prev2;

      return *this;
    }

    void SetJob(Job* &_job){
      job = _job;
      prev1.SetHeadNode(this);
      prev2.SetHeadNode(this);
    }

    int GetWeight(){
      return weight;
    }

    bool GetTerminalNode(){
      return terminal_node;
    }

    Job *GetJob(){
      return job;
    }

    friend bool operator<(const NodeDuration<T> &lhs, const NodeDuration<T> &rhs) {
        return lhs.prev1.GetF() < rhs.prev1.GetF();
    }

    friend bool operator> (const NodeDuration<T> lhs, const NodeDuration<T> rhs){ return rhs < lhs; }
    friend bool operator<=(const NodeDuration<T> lhs, const NodeDuration<T> rhs){ return !(lhs > rhs); }
    friend bool operator>=(const NodeDuration<T> lhs, const NodeDuration<T> rhs){ return !(lhs < rhs); }
};


#endif // NODE_DURATION_HPP
