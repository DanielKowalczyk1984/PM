#ifndef __LOCALSEARCH_NEW_H__
#define __LOCALSEARCH_NEW_H__

#include <array>
#include <list>
#include <vector>
#include "Solution_new.hpp"

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

    int operator()(int t) { return c + alpha * (t - b1); };
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
using SlopeTypeIt = std::list<SlopeType>::iterator;

struct LocalSearchData {

    LocalSearchData(int _nb_jobs, int _nb_machines);

    void RVND(Sol& sol);

    void swap_operator(Sol& sol, int l1, int l2);
    void swap_operator_inter(Sol& sol, int l1, int l2);
    void insertion_operator(Sol& sol, int l);
    void insertion_operator_inter(Sol& sol, int l);

   private:
    int  nb_jobs;
    int  nb_machines;
    bool updated{};

    MatrixInt                                      W;
    std::vector<std::vector<std::list<SlopeType>>> g;

    std::array<MatrixProcessingList, 2> processing_list;

    MatrixInt                G;
    MatrixInt                H;
    VectorInt                GG;
    MatrixInt                HH;
    std::vector<SlopeTypeIt> begin;
    std::vector<SlopeTypeIt> end;
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
    void calculate_g(const Sol& sol);
    void create_processing_insertion_operator(const Sol& sol, int l);
    void create_processing_insertion_operator_inter(const Sol& sol, int l);
    void create_processing_list_swap(const Sol& sol, int l1, int l2);
    void create_processing_list_swap_inter(const Sol& sol, int l1, int l2);

};

#endif  // __LOCALSEARCH_NEW_H__