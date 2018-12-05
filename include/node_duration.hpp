#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <memory>
#include <OptimalSolution.hpp>

template<typename T>
class Node;

template<typename T>
class PrevNode {
  private:
    PrevNode<T> *prev;
    bool high;
    Node<T>* head_node;

  public:
    T f;
    Job *prev_job;
    /**
     * Constructor
     */
    PrevNode(T &_f, PrevNode<T> *&_prev, bool &_high) :
      prev(_prev),
      high(_high),
      head_node(nullptr),
      f(_f),
      prev_job(nullptr){};

    PrevNode() :
      prev(nullptr),
      high(false),
      head_node(nullptr),
      f(-DBL_MAX),
      prev_job(nullptr){};

    /**
     * Copy Constructor
     */
    PrevNode<T>(const PrevNode<T> &src) :
      prev(src.prev),
      high(src.high),
      head_node(src.head_node),
      f(src.f),
      prev_job(src.prev_job) {}

    /**
     * Move Constructor
     */
    PrevNode<T>(const PrevNode<T> &&src) :
      prev(src.prev),
      high (src.high),
      head_node(src.head_node),
      f(src.f),
      prev_job(src.prev_job) {}

    /**
     * Copy Assignment
     */
    PrevNode<T>& operator=(const PrevNode<T> &src){
      if(&src == this) {
        return *this;
      }

      prev = src.prev;
      high = src.high;
      head_node = src.head_node;
      f = src.f;
      prev_job = src.prev_job;

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
      head_node = src.head_node;
      prev_job = src.prev_job;

      return *this;
    }

    void SetPrev(PrevNode<T> * &&_prev) {
        prev = _prev;
    }

    void SetF(T _f) {
        f = _f;
    }

    void SetHigh(bool &&_high) {
        high = _high;
    }

    void SetHeadNode(Node<T>* _head){
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

    bool GetHigh() {
        return high;
    }

    Job* GetJob() {
        return head_node->GetJob();
    }

    Job* GetPrevJob() {
      return GetPrev() == nullptr ? nullptr : GetPrev()->GetJob();
    }

    void UpdateSolution(T _f, PrevNode<T>* && _prev, bool &&_high) {
      f = _f;
      prev = _prev;
      high = _high;
    }

    void UpdateSolution(PrevNode<T> &_node) {
      f = _node.f;
      prev = _node.prev;
      high = _node.high;
    }

    Node<T> *GetNode() const {
      return head_node;
    }

    int GetWeight(){
      return head_node->GetWeight();
    }

    void UpdateNode(PrevNode<T> &_n){
      f = _n.f;
      high = _n.high;
      prev = _n.prev;
      prev_job = _n.prev_job;
    }

    void UpdateNode(T _f, Job *_job, bool &&_high) {
      f = _f;
      high = _high;
      prev_job = _job;
    }

    Job* get_prev_job() {
      return prev_job;
    }
};

template<typename T>
class Node {
  private:
    int weight;
    int num_layer;

    bool root_node;
    bool terminal_node;

    Job *job;

  public:
    PrevNode<T> prev1;
    PrevNode<T> prev2;

    std::shared_ptr<Node<T>> y;
    std::shared_ptr<Node<T>> n;

    Node<T>* child[2];
    
    /**
     * Constructor
     */
    Node():weight(0),
          num_layer(0),
          root_node(false),
          terminal_node(false),
          job(nullptr),
          prev1(),
          prev2(),
          y(nullptr),
          n(nullptr) {
        child[0] = nullptr;
        child[1] = nullptr;
        prev1.SetHeadNode(this);
        prev2.SetHeadNode(this);
    };

    Node(int &_weight, int &_num_layer, bool &_root_node,bool &_terminal_node):
         weight(_weight),
         num_layer(_num_layer),
         root_node(_root_node),
         terminal_node(_terminal_node),
         prev1(),
         prev2(),
         y(nullptr),
         n(nullptr) {
          prev1.SetHeadNode(this);
          prev2.SetHeadNode(this);
          child[0] = nullptr;
          child[1] = nullptr;
    }

    /**
     * Copy Constructor
     */
    Node<T>(const Node<T> &src):
        weight(src.weight),
        num_layer(src.num_layer),
        root_node(src.root_node),
        terminal_node(src.terminal_node),
        job(src.job),
        prev1(src.prev1),
        prev2(src.prev2),
        y(src.y),
        n(src.n),
        child{src.child[0], src.child[1]} {
      // child[0] = src.child[0];
      // child[1] = src.child[1]; 
      // child = {src.child[0], src.child[1]};
    }

    /**
     * Move Constructor
     */
    Node<T>(Node<T> &&src):
        weight(src.weight),
        num_layer(src.num_layer),
        root_node(src.root_node),
        terminal_node(src.terminal_node),
        job(src.job),
        prev1(src.prev1),
        prev2(src.prev2), 
        y(std::move(src.y)),
        n(std::move(src.n)){
        child = {src.child[0], src.child[1]};
      // child[0] = src.child[0];
      // child[1] = src.child[1];
    }

    /**
     * Copy Operator
     */
    Node<T>& operator=(const Node<T> &src){
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

      child = {src.child[0], src.child[1]};

      return *this;
    }

    /**
     * Move Operator
     */
    Node<T>& operator=(Node<T> &&src){
      if(&src == this) {
        return *this;
      }

      weight = src.weight;
      num_layer = src.num_layer;
      root_node = src.root_node;
      terminal_node = src.node;
      job = src.job;

      y = std::move(src.y);
      n = std::move(src.n);

      prev1 = src.prev1;
      prev2 = src.prev2;

      child = {src.child[0], src.child[1]};

      return *this;
    }

    // void SetJob(Job* &_job){
    //   job = _job;
    //   prev1.SetHeadNode(this);
    //   prev2.SetHeadNode(this);
    // }

    void set_job(Job *_job, bool _terminal_node = false){
      job = _job;
      terminal_node = _terminal_node;
    }

    void set_weight(int _weight){
      weight = _weight;
    }

    int GetWeight(){
      return weight;
    }

    int GetLayerNum() {
      return num_layer;
    }

    bool GetTerminalNode(){
      return terminal_node;
    }

    Job *GetJob(){
      return job;
    }

    Node<T>* InitNode(int _weight, bool _root_node = false, bool _terminal_node = false){
      if(!terminal_node) {
        weight = _weight;
        root_node = _root_node;
      } else {
        job = nullptr;
        weight = -1;
        root_node = _root_node;
        terminal_node = _terminal_node;
      }
      return this;
    }

    friend bool operator<(const Node<T> &lhs, const Node<T> &rhs) {
        return lhs.prev1.f < rhs.prev1.f;
    }

    friend bool operator> (const Node<T> &lhs, const Node<T> &rhs){ return rhs < lhs; }
    friend bool operator<=(const Node<T> &lhs, const Node<T> &rhs){ return !(lhs > rhs); }
    friend bool operator>=(const Node<T> &lhs, const Node<T> &rhs){ return !(lhs < rhs); }
};


#endif // NODE_DURATION_HPP
