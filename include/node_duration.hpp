#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <memory>
#include <OptimalSolution.hpp>

template<typename T>
class Node;

template<typename T>
class Label {
  private:
    Label<T> *prev;
    bool high;
    Node<T>* head_node;

  public:
    T f;
    Job *prev_job;
    /**
     * Constructor
     */
    Label(T &_f, Label<T> *&_prev, bool &_high) :
      prev(_prev),
      high(_high),
      head_node(nullptr),
      f(_f),
      prev_job(nullptr){};

    Label() :
      prev(nullptr),
      high(false),
      head_node(nullptr),
      f(-DBL_MAX),
      prev_job(nullptr){};

    /**
     * Copy Constructor
     */
    Label<T>(const Label<T> &src) :
      prev(src.prev),
      high(src.high),
      head_node(src.head_node),
      f(src.f),
      prev_job(src.prev_job) {}

    /**
     * Move Constructor
     */
    Label<T>(const Label<T> &&src) :
      prev(src.prev),
      high (src.high),
      head_node(src.head_node),
      f(src.f),
      prev_job(src.prev_job) {}

    /**
     * Copy Assignment
     */
    Label<T>& operator=(const Label<T> &src){
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
    Label<T>& operator=(const Label<T> &&src){
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

    void SetPrev(Label<T> * &&_prev) {
        prev = _prev;
    }

    void SetF(T _f) {
        f = _f;
    }

    void SetHigh(bool &&_high) {
        high = _high;
    }

    void SetHeadNode(Node<T>* _head) {
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

    Label<T>* GetPrev() {
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

    void UpdateSolution(T _f, Label<T>* && _prev, bool &&_high) {
      f = _f;
      prev = _prev;
      high = _high;
    }

    void UpdateSolution(Label<T> &_node) {
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

    void UpdateNode(Label<T> *_n, T _f = 0, bool _high = false){
      if(_high) {
        f = _f;
        prev_job = GetJob();
      } else {
        f = _n->f;
        prev_job = _n->prev_job;
      }
      high = _high;
      prev = _n;
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
    Label<T> forward_label1;
    Label<T> forward_label2;

    Label<T> backward_label1;
    Label<T> backward_label2;

    std::shared_ptr<Node<T>> y;
    std::shared_ptr<Node<T>> n;

    Node<T>* child[2];

    T dist_root_node;
    T dist_terminal_yes;
    T dist_terminal_no;

    bool calc_yes;
    bool calc_no;
    bool remove_node;
    
    /**
     * Constructor
     */
    Node():weight(0),
          num_layer(0),
          root_node(false),
          terminal_node(false),
          job(nullptr),
          forward_label1(),
          forward_label2(),
          backward_label1(),
          backward_label2(),
          y(nullptr),
          n(nullptr),
          calc_yes(true),
          calc_no(true),
          remove_node(false) {
        child[0] = nullptr;
        child[1] = nullptr;
        backward_label1.SetHeadNode(this);
        backward_label2.SetHeadNode(this);
        forward_label1.SetHeadNode(this);
        forward_label2.SetHeadNode(this);
    };

    Node(int &_weight, int &_num_layer, bool &_root_node,bool &_terminal_node):
         weight(_weight),
         num_layer(_num_layer),
         root_node(_root_node),
         terminal_node(_terminal_node),
         forward_label1(),
         forward_label2(),
         backward_label1(),
         backward_label2(),
         y(nullptr),
         n(nullptr),
         calc_yes(true),
         calc_no(true),
         remove_node(false) {
          forward_label1.SetHeadNode(this);
          forward_label2.SetHeadNode(this);
          backward_label1.SetHeadNode(this);
          backward_label2.SetHeadNode(this);
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
        forward_label1(src.forward_label1),
        forward_label2(src.forward_label2),
        backward_label1(src.backward_label1),
        backward_label2(src.backward_label2),
        y(src.y),
        n(src.n),
        child{src.child[0], src.child[1]},
        calc_yes(src.calc_yes),
        calc_no(src.calc_no),
        remove_node(src.remove_node) {
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
        forward_label1(src.forward_label1),
        forward_label2(src.forward_label2),
        backward_label1(src.backward_label1),
        backward_label2(src.backward_label2), 
        y(std::move(src.y)),
        n(std::move(src.n)),
        child{src.child[0], src.child[1]},
        calc_yes(src.calc_yes),
        calc_no(src.calc_no),
        remove_node(src.remove_node){
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

      forward_label1 = src.forward_label1;
      forward_label2 = src.forward_label2;

      backward_label1 = src.backward_label1;
      backward_label2 = src.backward_label2;

      child = {src.child[0], src.child[1]};

      dist_root_node = src.dist_root_node;
      dist_terminal_yes = src.dist_terminal_yes;
      dist_terminal_no = src.dist_terminal_no;

      calc_no = src.calc_no;
      calc_yes = src.calc_yes;
      remove_node = src.remove_node;

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

      forward_label1 = src.forward_label1;
      forward_label2 = src.forward_label2;

      backward_label1 = src.backward_label1;
      backward_label2 = src.backward_label2;

      child = {src.child[0], src.child[1]};
      dist_root_node = src.dist_root_node;
      dist_terminal_yes = src.dist_terminal_yes;
      dist_terminal_no = src.dist_terminal_no;

      calc_no = src.calc_no;
      calc_yes = src.calc_yes;
      remove_node = src.remove_node;
      return *this;
    }

    // void SetJob(Job* &_job){
    //   job = _job;
    //   state1.SetHeadNode(this);
    //   state2.SetHeadNode(this);
    // }

    void set_job(Job *_job, bool _terminal_node = false){
      job = _job;
      terminal_node = _terminal_node;
    }

    void set_weight(int _weight){
      weight = _weight;
    }

    void set_layer(int _num_layer) {
      num_layer = _num_layer;
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
      if(!_terminal_node) {
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
        return lhs.forward_label1.f < rhs.forward_label1.f;
    }

    friend bool operator> (const Node<T> &lhs, const Node<T> &rhs){ return rhs < lhs; }
    friend bool operator<=(const Node<T> &lhs, const Node<T> &rhs){ return !(lhs > rhs); }
    friend bool operator>=(const Node<T> &lhs, const Node<T> &rhs){ return !(lhs < rhs); }
};


#endif // NODE_DURATION_HPP
