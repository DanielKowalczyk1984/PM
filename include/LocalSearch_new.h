#ifndef __LOCALSEARCH_NEW_H__
#define __LOCALSEARCH_NEW_H__

#include <array>
#include <initializer_list>
#include <list>
#include <vector>
#include "Solution_new.hpp"

template <typename IT>
void swap_ranges(IT start_a, IT end_a, IT start_b, IT end_b) {
    // Will return the position after the moved elements
    auto it = std::rotate(start_a, start_b, end_b);
    // Will determine the point where the new range needs to be moved from
    auto new_start_a = (end_a - start_a) + it;
    std::rotate(it, new_start_a, end_b);
}

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
};

using MatrixProcessingList = std::vector<std::vector<ProcessingListData>>;

struct LocalSearchData {
    int  nb_jobs;
    int  nb_machines;
    bool updated{};

    std::vector<std::vector<int>>                  W;
    std::vector<std::vector<std::list<SlopeType>>> g;

    std::array<MatrixProcessingList, 2> processing_list;

    LocalSearchData(int _nb_jobs, int _nb_machines);

    void forward_insertion(Sol& sol, int l);
    void swap_operator(Sol& sol, int l1, int l2);
    void swap_operator_inter(Sol& sol, int l1, int l2);
    void insertion_operator_inter(Sol& sol, int l);

   private:
    void calculate_W(const Sol& sol);
    void calculate_g(const Sol& sol);
    void create_processing_list(const Sol& sol, int j);
    void create_processing_list_reversed(const Sol& sol, int j);
    void create_processing_list_swap(const Sol& sol, int l1, int l2);
    void create_processing_list_insertion(const Sol& sol, int l);
    void create_processing_list_swap_inter(const Sol& sol, int l1, int l2);

    void update_insertion(Sol& sol, int i_best, int j_best, int k_best, int l);
    void update_swap(Sol& sol,
                     int  i_best,
                     int  j_best,
                     int  k_best,
                     int  l1,
                     int  l2);
    void update_swap_inter(Sol& sol,
                           int  i_best,
                           int  j_best,
                           int  k_best,
                           int  kk_best,
                           int  l1,
                           int  l2);
    void update_insertion_inter(Sol& sol,
                                int  i_best,
                                int  j_best,
                                int  k_best,
                                int  kk_best,
                                int  l);

    std::vector<std::vector<int>>               G;
    std::vector<std::vector<int>>               H;
    std::vector<int>                            GG;
    std::vector<std::vector<int>>               HH;
    std::vector<std::list<SlopeType>::iterator> begin;
    std::vector<std::list<SlopeType>::iterator> end;
    std::vector<std::vector<int>>               B2_1;
    std::vector<std::vector<int>>               B2_2;
    std::vector<std::vector<int>>               B3_1;
    std::vector<std::vector<int>>               B3_2;
    std::vector<std::vector<int>>               B4_1;
    std::vector<std::vector<int>>               B4_2;
    std::vector<int>                            B3_1_;
    std::vector<std::vector<int>>               B5_1;
    std::vector<std::vector<int>>               B5_2;
    std::vector<std::vector<int>>               B6_1;
};

#endif  // __LOCALSEARCH_NEW_H__