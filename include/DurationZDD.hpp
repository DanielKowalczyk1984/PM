#ifndef DURATION_ZDD_HPP
#define DURATION_ZDD_HPP
#include <tdzdd/DdEval.hpp>

#include <OptimalSolution.hpp>
#include <node_duration.hpp>

#include <vector>
#include <algorithm>
using namespace std;

template<typename T>
class PricerDurationZDD {
  public:
    std::vector<std::shared_ptr<NodeDuration<T>>>  list;
    Job *job;

    PricerDurationZDD() {}

    ~PricerDurationZDD() {
        list.clear();
    }


    std::shared_ptr<NodeDuration<T>> add_weight(int _weight, int _layer, bool _rootnode = false, bool _terminal_node = false){
        for (auto &it : list) {
            if (it->GetWeight() == _weight) {
                return it;
            }
        }

        std::shared_ptr<NodeDuration<T>> node = std::make_shared<NodeDuration<T>>(_weight, _layer, _rootnode, _terminal_node);
        node->set_job(job);
        list.push_back(std::move(node));
        return list.back();
    }

    void set_job(Job *_job){
        job = _job;
    }

    Job* get_job(){
        return job;
    }
};

template<typename T>
bool my_compare(const std::shared_ptr<NodeDuration<T>> &lhs, const std::shared_ptr<NodeDuration<T>> &rhs){
    return *lhs < *rhs;
}

template<typename E, typename T> class DurationZDD : public
    tdzdd::DdEval<E, PricerDurationZDD<T>, Optimal_Solution<T>> {
  private:
    T *pi;
    int num_jobs;
  public:
    DurationZDD(T *_pi, int _num_jobs): pi(_pi), num_jobs(_num_jobs) {
    }

    DurationZDD(int _num_jobs)
    :  num_jobs(_num_jobs){
    }

    DurationZDD(){
        pi = nullptr;
        num_jobs = 0;
    }

    DurationZDD(const DurationZDD<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializenode(PricerDurationZDD<T>& n){
        for (auto &it: n.list) {
            if(it->GetWeight() == 0) {
                it->prev1.UpdateSolution(pi[num_jobs], nullptr, false);
                it->prev2.UpdateSolution(-DBL_MAX/2, nullptr, false);
            } else {
                it->prev1.UpdateSolution(-DBL_MAX/2, nullptr, false);
                it->prev2.UpdateSolution(-DBL_MAX/2, nullptr, false);
            }
        }
    }

    void initializerootnode(PricerDurationZDD<T> &n) {
        for(auto &it: n.list){
            it->prev1.f = pi[num_jobs];
            it->prev2.SetF(-DBL_MAX/2);
        }
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    void evalNode(PricerDurationZDD<T> &n) const
    {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);
        double result;

        for (auto &it : n.list) {
            int      weight = it->GetWeight();
            T g;
            std::shared_ptr<NodeDuration<T>> p0 = it->n;
            std::shared_ptr<NodeDuration<T>> p1 = it->y;
            result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

            /**
             * High edge calculation
             */
            Job *prev = it->prev1.GetPrevJob();
            Job *aux1 = p1->prev1.GetPrevJob();
            if(prev != tmp_j ) {
                g = it->prev1.GetF() + result;
                if(g > p1->prev1.GetF()) {
                    if(aux1 != tmp_j) {
                        p1->prev2.UpdateSolution(p1->prev1);
                    }
                    p1->prev1.UpdateSolution(g, &(it->prev1), true);
                } else if ((g > p1->prev2.GetF()) && (aux1 != tmp_j)) {
                    p1->prev2.UpdateSolution(g, &(it->prev1), true);
                }
            } else  {
                g = it->prev2.GetF() + result;
                if(g > p1->prev1.GetF()) {
                    if(aux1 != tmp_j) {
                        p1->prev2.UpdateSolution(p1->prev1);
                    }
                    p1->prev1.UpdateSolution(g, &(it->prev2), true);
                } else if ((g >= p1->prev2.GetF()) && (aux1 != tmp_j)) {
                    p1->prev2.UpdateSolution(g, &(it->prev2), true);
                }
            }

            /**
             * Low edge calculation
             */
            aux1 = p0->prev1.GetPrevJob();
            if(it->prev1.GetF() > p0->prev1.GetF()) {
                if(prev != aux1) {
                    p0->prev2.UpdateSolution(p0->prev1);
                }
                p0->prev1.UpdateSolution(it->prev1);
            } else if ((it->prev1.GetF() > p0->prev2.GetF()) && (aux1 != prev)){
                p0->prev2.UpdateSolution(it->prev1);
            }
        }
    }


    Optimal_Solution<T> get_objective(PricerDurationZDD<T> &n) {
        Optimal_Solution<T> sol;
        Job *aux_job;

        sol.cost = 0;
        sol.C_max = 0;
        sol.obj = pi[num_jobs];
        int weight;


        auto m =  std::max_element(n.list.begin(), n.list.end(),my_compare<T>);
        weight = (*m)->GetWeight();

        PrevNode<T> *ptr_node = &((*m)->prev1);

        while(ptr_node->GetPrev() != nullptr) {
            PrevNode<T> *aux_prev_node = ptr_node->GetPrev();
            aux_job = aux_prev_node->GetJob();
            if(ptr_node->GetHigh()) {
                sol.C_max += aux_job->processingime;
                if(ptr_node->GetNode()) {
                    sol.cost = sol.cost + value_Fj(ptr_node->GetWeight(), aux_job) ;
                    sol.obj += pi[aux_job->job] - value_Fj(ptr_node->GetWeight(), aux_job);
                } else {
                    sol.cost = sol.cost + value_Fj(weight, aux_job);
                    sol.obj += pi[aux_job->job] - value_Fj(weight, aux_job);
                }
                g_ptr_array_insert(sol.jobs, 0, aux_job);
                // g_ptr_array_add(sol.jobs, aux_job);
            }
            ptr_node = aux_prev_node;
        }

        assert(sol.C_max == weight);
        // printf("test cost %f %f\n", sol.obj, (*m)->prev1.GetF());
        // for(int i = 0; i < sol.jobs->len; i++){
        //     Job *tmp_job = (Job *) g_ptr_array_index(sol.jobs, i);
        //     printf("%d ", tmp_job->job);
        // }
        // printf("\n");

        return sol;
    }
};

#endif // DURATION_ZDD_HPP