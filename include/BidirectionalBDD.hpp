#ifndef BIDIRECTIONAL_BDD_HPP
#define BIDIRECTIONAL_BDD_HPP
#include <BackwardBDD.hpp>
#include <ForwardBDD.hpp>

template<typename E, typename T> class BidirectionalBddSimple : public 
    tdzdd::DdEval<E, Node<T>, Optimal_Solution<T>> {

    struct ForwardBddSimpleBidirectional : ForwardBddSimple<ForwardBddSimpleBidirectional, T> {
        ForwardBddSimpleBidirectional(T *_pi, int _num_jobs)
        : ForwardBddSimple<ForwardBddSimpleBidirectional, T>(_pi, _num_jobs) {};
        explicit ForwardBddSimpleBidirectional(int _num_jobs)
        : ForwardBddSimple<ForwardBddSimpleBidirectional, T> (_num_jobs) {};
        ForwardBddSimpleBidirectional() : ForwardBddSimple<ForwardBddSimpleBidirectional, T>() {};
    };

    struct BackwardBddSimpleBidirectional : BackwardBddSimple<BackwardBddSimpleBidirectional, T> {
        BackwardBddSimpleBidirectional(T *_pi, int _num_jobs)
        : BackwardBddSimple<BackwardBddSimpleBidirectional, T>(_pi, _num_jobs) {};
        explicit BackwardBddSimpleBidirectional(int _num_jobs)
        : BackwardBddSimple<BackwardBddSimpleBidirectional, T> (_num_jobs) {};
        BackwardBddSimpleBidirectional() : BackwardBddSimple<BackwardBddSimpleBidirectional, T>() {};
    };

private:
    ForwardBddSimpleBidirectional forward_evaluator;
    BackwardBddSimpleBidirectional backward_evaluator;
    int num_layer;
public:
    BidirectionalBddSimple(T *_pi, int _num_jobs, int _num_layer): forward_evaluator(_pi,_num_jobs), backward_evaluator(_pi,_num_jobs), num_layer(_num_layer) {
    }

    BidirectionalBddSimple(int _num_jobs, int _num_layer)
    :  forward_evaluator(_num_jobs), backward_evaluator(_num_jobs), num_layer(_num_layer){
    
    }

    BidirectionalBddSimple() : forward_evaluator(), backward_evaluator(), num_layer{}{
    }

    BidirectionalBddSimple(const BidirectionalBddSimple<E, T> &src) {
        forward_evaluator(src.forward_evaluator);
        backward_evaluator(src.backward_evaluator);
        num_layer = src.num_layer;
    }

    void initializepi(T *_pi){
        forward_evaluator.initializepi(_pi);
        backward_evaluator.initializepi(_pi);
    }

    void initializenode(Node<T>& n) {
        if (n.GetLayerNum() < num_layer) {
            forward_evaluator.initializenode(n);
        } else {
            backward_evaluator.initializenode(n);
        }

    };

    void initializerootnode(Node<T>& n) const {
        forward_evaluator.initializerootnode(n);
    };

    void evalNode(Node<T>& n) const {
        if(n.GetLayerNum() < num_layer) {
            forward_evaluator.evalNode(n);
        } else {
            backward_evaluator.evalNode(n);
        }
    };

    Optimal_Solution<T> get_objective(Node<T> &n) {
        Optimal_Solution<T> sol;

        return sol;
    };
};


#endif // BIDIRECTIONAL_BDD_HPP
