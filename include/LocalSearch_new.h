#ifndef __LOCALSEARCH_NEW_H__
#define __LOCALSEARCH_NEW_H__

#include <array>
#include <boost/timer/timer.hpp>
#include <list>
#include <random>
#include <vector>
#include "Solution_new.hpp"
#include "util.h"

struct SlopeType {
    int b1{};    /* begin of segment*/
    int b2{};    /* end of segment*/
    int c{};     /* value of the function at b1*/
    int alpha{}; /* slope of the function */

    SlopeType() = default;
    SlopeType(const SlopeType&) = default;
    SlopeType(SlopeType&&) = default;
    SlopeType& operator=(SlopeType&&) = default;
    SlopeType& operator=(const SlopeType&) = default;
    ~SlopeType() = default;

    SlopeType(int _b1, int _b2, int _c, int _alpha)
        : b1(_b1),
          b2(_b2),
          c(_c),
          alpha(_alpha) {}

    int  operator()(int t) { return c + alpha * (t - b1); };
    bool in_interval(int c) const { return (b1 <= c) && (c < b2); };
};

struct ProcessingListData {
    std::size_t pos{};
    int         p{};

    ProcessingListData() = default;
    ProcessingListData(std::size_t _pos, int _p) : pos(_pos), p(_p){};
    ProcessingListData& operator=(const ProcessingListData&) = default;
    ProcessingListData(const ProcessingListData&) = default;
    ProcessingListData(ProcessingListData&&) = default;
    ProcessingListData& operator=(ProcessingListData&&) = default;
    ~ProcessingListData() = default;
};

using MatrixProcessingList = std::vector<std::vector<ProcessingListData>>;
using VectorProcessingList = std::vector<ProcessingListData>;
using MatrixInt = std::vector<std::vector<int>>;
using VectorInt = std::vector<int>;
using SlopeList = std::list<SlopeType>;
using SlopeListIt = std::list<SlopeType>::iterator;

struct LocalSearchData {
    LocalSearchData(int _nb_jobs, int _nb_machines);

    void RVND(Sol& sol);

    void swap_operator(Sol& sol, int l1, int l2);
    void swap_operator_inter(Sol& sol, int l1, int l2);
    void insertion_operator(Sol& sol, int l);
    void insertion_operator_inter(Sol& sol, int l);
    int  get_iterations() const { return iterations; }

    boost::timer::cpu_timer test_swap;
    CCutil_timer            test_timer;

   private:
    int  nb_jobs;
    int  nb_machines;
    bool updated{};
    int  iterations{};

    MatrixInt                           W;
    std::vector<std::vector<SlopeList>> g;

    std::array<MatrixProcessingList, 2> processing_list;

    MatrixInt                G;
    MatrixInt                H;
    VectorInt                GG;
    MatrixInt                HH;
    std::vector<SlopeListIt> begin;
    std::vector<SlopeListIt> end;
    MatrixInt                B2_1;
    MatrixInt                B2_2;
    MatrixInt                B3_1;
    MatrixInt                B3_2;
    MatrixInt                B4_1;
    MatrixInt                B4_2;
    VectorInt                B3_1_;
    MatrixInt                B5_1;
    MatrixInt                B5_2;
    MatrixInt                B6_1;

    void calculate_W(const Sol& sol);
    void calculate_g(Sol& sol);
    void create_processing_insertion_operator(const Sol& sol, int l);
    void create_processing_insertion_operator_inter(const Sol& sol, int l);
    void create_processing_list_swap(const Sol& sol, int l1, int l2);
    void create_processing_list_swap_inter(const Sol& sol, int l1, int l2);
};

struct MoveOperator {
    LocalSearchData* data{};

    explicit MoveOperator(LocalSearchData* _data) : data(_data){};
    MoveOperator() = default;
    MoveOperator(const MoveOperator&) = default;
    MoveOperator(MoveOperator&&) = default;
    MoveOperator& operator=(MoveOperator&&) = default;
    MoveOperator& operator=(const MoveOperator&) = default;
    virtual ~MoveOperator() = default;
    virtual bool operator()(Sol& sol) = 0;

    int get_iterations() const { return data->get_iterations(); }
};

struct SwapOperator : public MoveOperator {
    int l1{};
    int l2{};

    SwapOperator(LocalSearchData* _data, int _l1, int _l2)
        : MoveOperator(_data),
          l1(_l1),
          l2(_l2){};
    SwapOperator() = default;
    SwapOperator(const SwapOperator&) = default;
    SwapOperator(SwapOperator&&) = default;
    SwapOperator& operator=(SwapOperator&&) = default;
    SwapOperator& operator=(const SwapOperator&) = default;
    virtual ~SwapOperator() = default;

    bool operator()(Sol& sol) override {
        data->swap_operator(sol, l1, l2);
        return true;
    }
};

struct SwapOperatorInter : public MoveOperator {
    int l1{};
    int l2{};

    SwapOperatorInter(LocalSearchData* _data, int _l1, int _l2)
        : MoveOperator(_data),
          l1(_l1),
          l2(_l2){};
    SwapOperatorInter() = default;
    SwapOperatorInter(const SwapOperatorInter&) = default;
    SwapOperatorInter(SwapOperatorInter&&) = default;
    SwapOperatorInter& operator=(SwapOperatorInter&&) = default;
    SwapOperatorInter& operator=(const SwapOperatorInter&) = default;
    virtual ~SwapOperatorInter() = default;

    bool operator()(Sol& sol) override {
        data->swap_operator_inter(sol, l1, l2);
        return true;
    }
};

struct InsertionOperator : public MoveOperator {
    int l{};

    InsertionOperator(LocalSearchData* _data, int _l)
        : MoveOperator(_data),
          l(_l){};
    InsertionOperator() = default;
    InsertionOperator(const InsertionOperator&) = default;
    InsertionOperator(InsertionOperator&&) = default;
    InsertionOperator& operator=(InsertionOperator&&) = default;
    InsertionOperator& operator=(const InsertionOperator&) = default;
    virtual ~InsertionOperator() = default;

    bool operator()(Sol& sol) override {
        data->insertion_operator(sol, l);
        return true;
    }
};

struct InsertionOperatorInter : public MoveOperator {
    int l{};

    InsertionOperatorInter(LocalSearchData* _data, int _l)
        : MoveOperator(_data),
          l(_l){};
    InsertionOperatorInter() = default;
    InsertionOperatorInter(const InsertionOperatorInter&) = default;
    InsertionOperatorInter(InsertionOperatorInter&&) = default;
    InsertionOperatorInter& operator=(InsertionOperatorInter&&) = default;
    InsertionOperatorInter& operator=(const InsertionOperatorInter&) = default;
    virtual ~InsertionOperatorInter() = default;

    bool operator()(Sol& sol) override {
        data->insertion_operator_inter(sol, l);
        return true;
    }
};

static std::mt19937 gen = std::mt19937(0);

struct PerturbationOperator {
    int l1{};
    int l2{};

    PerturbationOperator() = default;
    PerturbationOperator(int _l1, int _l2) : l1(_l1), l2(_l2) {}
    PerturbationOperator(const PerturbationOperator&) = default;
    PerturbationOperator(PerturbationOperator&&) = default;
    PerturbationOperator& operator=(const PerturbationOperator&) = default;
    PerturbationOperator& operator=(PerturbationOperator&&) = default;
    ~PerturbationOperator() = default;

    void operator()(Sol& sol) const { sol.perturb_solution(l1, l2, gen); }
};

#endif  // __LOCALSEARCH_NEW_H__