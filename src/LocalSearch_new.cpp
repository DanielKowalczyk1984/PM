#include "LocalSearch_new.h"
#include <bits/c++config.h>
#include <fmt/core.h>
#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include "Solution_new.hpp"
#include "job.h"

LocalSearchData::LocalSearchData(int _nb_jobs, int _nb_machines)
    : nb_jobs(_nb_jobs),
      nb_machines(_nb_machines),
      W(_nb_machines, VectorInt(_nb_jobs)),
      g(_nb_machines,
        std::vector<std::list<SlopeType>>(nb_jobs, std::list<SlopeType>())),
      processing_list{
          MatrixProcessingList(
              _nb_machines,
              VectorProcessingList(_nb_jobs, ProcessingListData())),
          MatrixProcessingList(
              _nb_machines,
              VectorProcessingList(_nb_jobs, ProcessingListData()))},
      G(_nb_jobs, VectorInt(_nb_jobs, 0)),
      H(_nb_jobs, VectorInt(_nb_jobs, 0)),
      GG(_nb_jobs, 0),
      HH(_nb_jobs, std::vector<int>(_nb_jobs, 0)),
      begin(_nb_jobs, std::list<SlopeType>::iterator()),
      end(_nb_jobs, std::list<SlopeType>::iterator()),
      B2_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B2_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B4_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B4_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_1_(_nb_jobs + 1, 0),
      B5_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B5_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B6_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)) {}

void LocalSearchData::RVND(Sol& sol) {
    do {
    } while (updated);
}

void LocalSearchData::insertion_operator(Sol& sol, int l) {
    int max = 0;
    updated = false;
    auto i_best{0};
    auto j_best{0};
    auto k_best{0};

    calculate_g(sol);
    calculate_W(sol);
    create_processing_insertion_operator(sol, l);

    for (auto k = 0UL; auto& m : sol.machines) {
        const auto& vec_m = m.job_list;
        const auto& n = vec_m.size();

        for (auto i = 0UL; i < n - l; i++) {
            auto&      p_list = processing_list[0][k][i];
            auto       it = g[k][p_list.pos].begin();
            const auto end = g[k][p_list.pos].end();

            for (auto j = p_list.pos + l; j < n; j++) {
                auto& tmp_j = vec_m[j];
                auto  c = sol.c[tmp_j->job] - p_list.p;
                it = std::find_if(it, end, [&c](const SlopeType& tmp) -> bool {
                    return (tmp.b1 <= c && tmp.b2 > c);
                });
                G[p_list.pos][j] = it != end ? (*it)(c) : 0;
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            auto       it = g[k][i + l].begin();
            const auto end = g[k][i + l].end();

            for (auto j = i + l; j < n; j++) {
                auto& tmp_j = vec_m[j];
                auto  c = sol.c[tmp_j->job];
                it = std::find_if(it, end, [&c](const SlopeType& tmp) -> bool {
                    return (tmp.b1 <= c && tmp.b2 > c);
                });
                H[i][j] = it != end ? (*it)(c) : 0;
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            auto       it = g[k][i + l].begin();
            const auto end = g[k][i + l].end();
            auto       c = 0;

            if (i != 0) {
                auto tmp_j = vec_m[i - 1];
                c = sol.c[tmp_j->job];
            }
            it = std::find_if(it, end, [&c](const SlopeType& tmp) -> bool {
                return (tmp.b1 <= c && tmp.b2 > c);
            });
            GG[i] = it != end ? (*it)(c) : 0;
        }

        for (auto j = l; j < n - 1; ++j) {
            begin[j] = g[k][j + 1].begin();
            end[j] = g[k][j + 1].end();
        }

        for (int i = 0; i < n - l; ++i) {
            auto pos = processing_list[0][k][i].pos;
            auto p = processing_list[0][k][i].p;

            for (int j = pos + l; j < n; ++j) {
                if (j != n - 1) {
                    auto& tmp_j = vec_m[j];
                    auto  c = sol.c[tmp_j->job] - p;
                    begin[j] = std::find_if(
                        begin[j], end[j], [&c](const SlopeType& tmp) -> bool {
                            return (tmp.b1 <= c && tmp.b2 > c);
                        });
                    HH[pos][j] = begin[j] != end[j] ? (*begin[j])(c) : 0;
                } else {
                    HH[pos][j] = 0;
                }
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            for (auto j = i + l; j < n; j++) {
                auto tmp = 0;

                if (i != 0) {
                    tmp += W[k][i - 1];
                }

                tmp += G[i][j] - H[i][j];
                tmp += GG[i] - HH[i][j];
                tmp += W[k][n - 1] - W[k][j];

                if (m.total_weighted_tardiness - tmp > max) {
                    max = m.total_weighted_tardiness - tmp;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    updated = true;
                }
            }
        }

        k++;
    }

    if (updated) {
        sol.update_insertion_move(i_best, j_best, k_best, l);
    }
}

void LocalSearchData::swap_operator(Sol& sol, int l1, int l2) {
    std::shared_ptr<Job> tmp_j;

    auto c = 0, t = 0;
    auto max = 0;
    auto i_best = -1, j_best = -1, k_best = -1;

    updated = false;
    calculate_g(sol);
    calculate_W(sol);
    create_processing_list_swap(sol, l1, l2);

    for (auto k = 0UL; auto& m : sol.machines) {
        auto&       machine = m.job_list;
        const auto& n = machine.size();

        /** compute g */
        for (auto i = 0UL; i < n - l1 - l2 + 1; ++i) {
            auto&      p_list = processing_list[0][k][i];
            auto       it = g[k][p_list.pos].begin();
            const auto end = g[k][p_list.pos].end();

            for (auto j = p_list.pos + l1; j < n - l2 + 1; ++j) {
                tmp_j = machine[j + l2 - 1];
                c = sol.c[tmp_j->job] - p_list.p;
                it = std::find_if(it, end, [&c](const auto& tmp) -> bool {
                    return (tmp.b1 <= c) && (tmp.b2 > c);
                });
                B2_1[p_list.pos][j] = it != end ? (*it)(c) : 0;
            }
        }

        /** compute h */
        for (int i = 0; i < n - l1 - l2 + 1; ++i) {
            if (i + l1 >= n) {
                for (int j = i + l1; j < n - l2 + 1; ++j) {
                    B2_2[i][j] = 0;
                }
            } else {
                auto       it = g[k][i + l1].begin();
                const auto end = g[k][i + l1].end();
                for (int j = i + l1; j < n - l2 + 1; ++j) {
                    tmp_j = machine[j + l2 - 1];
                    c = sol.c[tmp_j->job];
                    it = std::find_if(it, end, [&c](const auto& tmp) -> bool {
                        return (tmp.b1 <= c) && (tmp.b2 > c);
                    });
                    B2_2[i][j] = it != end ? (*it)(c) : 0;
                }
            }
        }

        /** compute B3_1 */
        for (int i = 0; i < n - l1 - l2 + 1; ++i) {
            begin[i] = g[k][i + l1].begin();
            end[i] = g[k][i + l1].end();
        }

        for (int j = l1; j < n - l2 + 1; ++j) {
            auto& p_list = processing_list[1][k][j - l1];

            for (int i = 0; i < p_list.pos - l1 + 1; ++i) {
                if (i + l1 >= n) {
                    B3_1[i][p_list.pos] = 0;
                } else {
                    c = p_list.p;

                    if (i != 0) {
                        tmp_j = machine[i - 1];
                        c += sol.c[tmp_j->job];
                    }

                    begin[i] = std::find_if(
                        begin[i], end[i], [&c](const auto& tmp) -> bool {
                            return (tmp.b1 <= c) && (tmp.b2 > c);
                        });
                    B3_1[i][p_list.pos] =
                        begin[i] != end[i] ? (*begin[i])(c) : 0;
                }
            }
        }

        /** compute B3_2 */
        for (int j = l1; j < n - l2 + 1; ++j) {
            begin[j] = g[k][j].begin();
            end[j] = g[k][j].end();
        }

        for (int i = 0; i < n - l1 - l2 + 1; ++i) {
            auto& p_list = processing_list[0][k][i];

            for (int j = p_list.pos + l1; j < n - l2 + 1; ++j) {
                if (i + l1 >= n) {
                    B3_2[p_list.pos][j] = 0;
                } else {
                    tmp_j = machine[j + l2 - 1];
                    c = sol.c[tmp_j->job] - p_list.p;
                    begin[j] = std::find_if(
                        begin[j], end[j], [&c](const auto& tmp) -> bool {
                            return (tmp.b1 <= c) && (tmp.b2 > c);
                        });
                    B3_2[p_list.pos][j] =
                        begin[j] != end[j] ? (*begin[j])(c) : 0;
                }
            }
        }

        /** compute B4_1 */
        for (int j = l1; j < n - l2 + 1; ++j) {
            auto       it = g[k][j].begin();
            const auto end = g[k][j].end();

            for (int i = 0; i < j - l1 + 1; ++i) {
                c = 0;

                if (i != 0) {
                    tmp_j = machine[i - 1];
                    c = sol.c[tmp_j->job];
                }

                it = std::find_if(it, end, [&c](const auto& tmp) -> bool {
                    return (tmp.b1 <= c) && (tmp.b2 > c);
                });
                B4_1[i][j] = it != end ? (*it)(c) : 0;
            }
        }

        /** compute B4_2 */
        for (int j = l1; j < n - l2 + 1; ++j) {
            auto& p_list = processing_list[1][k][j - l1];

            if (p_list.pos + l2 == n) {
                for (int i = 0; i < p_list.pos - l1 + 1; ++i) {
                    B4_2[i][p_list.pos] = 0;
                }
            } else {
                auto       it = g[k][p_list.pos + l2].begin();
                const auto end = g[k][p_list.pos + l2].end();

                for (int i = 0; i < p_list.pos - l1 + 1; ++i) {
                    c = p_list.p;

                    if (i != 0) {
                        tmp_j = machine[i - 1];
                        c += sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) -> bool {
                        return (tmp.b1 <= c) && (tmp.b2 > c);
                    });

                    B4_2[i][p_list.pos] = it != end ? (*it)(c) : 0;
                }
            }
        }

        for (int i = 0; i < n - l1 - l2 + 1; ++i) {
            for (int j = i + l1; j < n - l2 + 1; ++j) {
                t = 0;

                if (i != 0) {
                    t += W[k][i - 1];
                }

                t += B2_1[i][j] - B2_2[i][j];
                t += B3_1[i][j] - B3_2[i][j];
                t += B4_1[i][j] - B4_2[i][j];
                t += W[k][n - 1] - W[k][j + l2 - 1];

                if (m.total_weighted_tardiness - t > max) {
                    max = m.total_weighted_tardiness - t;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    updated = true;
                }
            }
        }
        k++;
    }

    /** update to best improvement */
    if (updated) {
        if (dbg_lvl()) {
            sol.print_solution();
        }

        sol.update_swap_move(i_best, j_best, k_best, l1, l2);

        if (dbg_lvl()) {
            sol.print_solution();
        }
    }

    if (dbg_lvl() > 0) {
        printf(
            "intra swap with l1 = %d and l2 = %d, running time = %f and "
            "improvement %d on machine %d on places %d %d\n",
            l1, l2, CCutil_zeit(), max, k_best, i_best, j_best);
    }
}

void LocalSearchData::swap_operator_inter(Sol& sol, int l1, int l2) {
    std::shared_ptr<Job> tmp_j;

    int    c, t;
    auto   update = false;
    auto   max = 0;
    auto   i_best = -1, j_best = -1, k_best = -1, kk_best = -1;
    double runningtime = CCutil_zeit();

    calculate_g(sol);
    calculate_W(sol);
    create_processing_list_swap_inter(sol, l1, l2);

    for (auto k1 = 0UL; k1 < nb_machines; ++k1) {
        auto&      machine1 = sol.machines[k1].job_list;
        const auto nb_jobs1 = machine1.size();

        for (auto k2 = 0UL; k2 < nb_machines; ++k2) {
            auto&      machine2 = sol.machines[k2].job_list;
            const auto nb_jobs2 = machine2.size();

            if (k1 == k2 || nb_jobs1 - l1 < 0 || nb_jobs2 - l2 < 0) {
                continue;
            }

            /** compute B2_1 */
            for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                auto       it = g[k1][i].begin();
                const auto end = g[k1][i].end();

                for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                    c = 0;

                    if (j != 0) {
                        tmp_j = machine2[j - 1];
                        c = sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) {
                        return (tmp.b1 <= c) && (c < tmp.b2);
                    });
                    B2_1[i][j] = it != end ? (*it)(c) : 0;
                }
            }

            for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                if (p_list.pos + l1 >= nb_jobs1) {
                    for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                        B2_2[p_list.pos][j] = 0;
                    }
                } else {
                    auto       it = g[k1][p_list.pos + l1].begin();
                    const auto end = g[k1][p_list.pos + l1].end();
                    for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                        c = p_list.p;

                        if (j != 0) {
                            tmp_j = machine2[j - 1];
                            c += sol.c[tmp_j->job];
                        }

                        it = std::find_if(it, end, [&c](const auto& tmp) {
                            return (tmp.b1 <= c) && (c < tmp.b2);
                        });
                        B2_2[p_list.pos][j] = it != end ? (*it)(c) : 0;
                    }
                }
            }

            for (int i = 0; i < nb_jobs1 - l1; ++i) {
                begin[i] = g[k1][i + l1].begin();
                end[i] = g[k1][i + l1].end();
            }

            for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                auto& p_list = processing_list[1][k2][j];

                for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                    if (i + l1 >= nb_jobs1) {
                        B3_1[i][p_list.pos] = 0;
                    } else {
                        c = p_list.p;

                        if (i != 0) {
                            tmp_j = machine1[i - 1];
                            c += sol.c[tmp_j->job];
                        }

                        begin[i] = std::find_if(
                            begin[i], end[i], [&c](const auto& tmp) {
                                return (tmp.b1 <= c) && (c < tmp.b2);
                            });
                        B3_1[i][p_list.pos] =
                            begin[i] != end[i] ? (*begin[i])(c) : 0;
                    }
                }
            }

            for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                auto it = g[k2][j].begin();
                auto end = g[k2][j].end();

                for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                    c = 0;

                    if (i != 0) {
                        tmp_j = machine1[i - 1];
                        c += sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) {
                        return (tmp.b1 <= c) && (c < tmp.b2);
                    });
                    B5_1[i][j] = it != end ? (*it)(c) : 0;
                }
            }

            for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                auto& p_list = processing_list[1][k2][j];

                if (p_list.pos + l2 >= nb_jobs2) {
                    for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                        B5_2[i][p_list.pos] = 0;
                    }
                } else {
                    auto it = g[k2][p_list.pos + l2].begin();
                    auto end = g[k2][p_list.pos + l2].end();
                    for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                        c = p_list.p;

                        if (i != 0) {
                            tmp_j = machine1[i - 1];
                            c += sol.c[tmp_j->job];
                        }

                        it = std::find_if(it, end, [&c](const auto& tmp) {
                            return (tmp.b1 <= c) && (c < tmp.b2);
                        });
                        B5_2[i][p_list.pos] = it != end ? (*it)(c) : 0;
                    }
                }
            }

            for (int j = 0; j < nb_jobs2 - l2; ++j) {
                begin[j] = g[k2][j + l2].begin();
                end[j] = g[k2][j + l2].end();
            }

            for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                    if (j + l2 >= nb_jobs2) {
                        B6_1[p_list.pos][j] = 0;
                    } else {
                        c = p_list.p;

                        if (j != 0) {
                            tmp_j = machine2[j - 1];
                            c += sol.c[tmp_j->job];
                        }

                        begin[j] = std::find_if(
                            begin[j], end[j], [&c](const auto& tmp) {
                                return (tmp.b1 <= c) && (tmp.b2 > c);
                            });
                        B6_1[p_list.pos][j] =
                            begin[j] != end[j] ? (*begin[j])(c) : 0;
                    }
                }
            }

            for (int i = 0; i < nb_jobs1 - l1 + 1; ++i) {
                for (int j = 0; j < nb_jobs2 - l2 + 1; ++j) {
                    t = 0;

                    if (i != 0) {
                        t += W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1[i][j];
                    t += B5_1[i][j] - B5_2[i][j];
                    t += B6_1[i][j];

                    if (j != 0) {
                        t += W[k2][j - 1];
                    }

                    if (sol.machines[k1].total_weighted_tardiness +
                            sol.machines[k2].total_weighted_tardiness - t >
                        max) {
                        max = sol.machines[k1].total_weighted_tardiness +
                              sol.machines[k2].total_weighted_tardiness - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        update = true;
                    }
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl()) {
            sol.print_solution();
        }

        sol.update_swap_move_inter(i_best, j_best, k_best, kk_best, l1, l2);
        updated = true;

        if (dbg_lvl()) {
            sol.print_solution();
        }
    } else {
        updated = false;
    }

    if (dbg_lvl() > 0) {
        printf(
            "inter insertion with l1 = %d and l2 = %d, running time = %f and "
            "improvement %d on machines %d and %d on places %d %d\n",
            l1, l2, CCutil_zeit() - runningtime, max, k_best, kk_best, i_best,
            j_best);
        print_line();
    }
}

void LocalSearchData::insertion_operator_inter(Sol& sol, int l) {
    int                  c, t;
    int                  update;
    std::shared_ptr<Job> tmp_j;
    int                  max;
    int                  i_best = -1, j_best = -1, k_best = -1, kk_best = -1;
    double               runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    calculate_g(sol);
    calculate_W(sol);
    create_processing_insertion_operator_inter(sol, l);

    for (auto k1 = 0UL; auto& m1 : sol.machines) {
        auto&       vec_m1 = m1.job_list;
        const auto& nb_jobs1 = vec_m1.size();

        for (auto k2 = 0UL; auto& m2 : sol.machines) {
            if (k1 == k2) {
                k2++;
                continue;
            }

            auto&       vec_m2 = m2.job_list;
            const auto& nb_jobs2 = vec_m2.size();

            for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                auto       it = g[k1][i].begin();
                const auto end = g[k1][i].end();

                for (auto j = 0UL; j < nb_jobs2; ++j) {
                    c = 0;

                    if (j != 0) {
                        tmp_j = vec_m2[j - 1];
                        c = sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) {
                        return (tmp.b1 <= c) && (tmp.b2 > c);
                    });
                    B2_1[i][j] = (it != end) ? (*it)(c) : 0;
                }
            }

            for (int i = 0; i < nb_jobs1 - l + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                if (p_list.pos + l >= nb_jobs1) {
                    for (int j = 0; j < nb_jobs2; ++j) {
                        B2_2[p_list.pos][j] = 0;
                    }
                } else {
                    auto       it = g[k1][p_list.pos + l].begin();
                    const auto end = g[k1][p_list.pos + l].end();
                    for (int j = 0UL; j < nb_jobs2; ++j) {
                        c = p_list.p;

                        if (j != 0UL) {
                            tmp_j = vec_m2[j - 1];
                            c += sol.c[tmp_j->job];
                        }

                        it = std::find_if(it, end, [&c](const auto& tmp) {
                            return (tmp.b1 <= c) && (tmp.b2 > c);
                        });
                        B2_2[p_list.pos][j] = it != end ? (*it)(c) : 0;
                    }
                }
            }

            for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                if (i + l >= nb_jobs1) {
                    B3_1_[i] = 0;
                } else {
                    auto       it = g[k1][i].begin();
                    const auto end = g[k1][i].end();
                    c = 0;

                    if (i != 0) {
                        tmp_j = vec_m1[i - 1];
                        c = sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) {
                        return (tmp.b1 <= c) && (tmp.b2 > c);
                    });
                    B3_1_[i] = it != end ? (*it)(c) : 0;
                }
            }

            for (int j = 0; j < nb_jobs2 - 1; ++j) {
                auto       it = g[k2][j].begin();
                const auto end = g[k2][j].end();

                for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                    auto& p_list = processing_list[0][k1][i];
                    c = p_list.p;

                    if (j != 0) {
                        tmp_j = vec_m2[j - 1];
                        c += sol.c[tmp_j->job];
                    }

                    it = std::find_if(it, end, [&c](const auto& tmp) {
                        return (tmp.b1 <= c) && (c < tmp.b2);
                    });
                    B5_1[p_list.pos][j] = (it != end) ? (*it)(c) : 0;
                }
            }

            for (int i = 0; i < nb_jobs1 - l + 1; ++i) {
                for (int j = 0; j < nb_jobs2 - 1; ++j) {
                    t = 0;

                    if (i != 0) {
                        t += W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1_[i];
                    t += B5_1[i][j];

                    if (j != 0) {
                        t += W[k2][j - 1];
                    }

                    if (m1.total_weighted_tardiness +
                            m2.total_weighted_tardiness - t >
                        max) {
                        max = m1.total_weighted_tardiness +
                              m2.total_weighted_tardiness - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        update = 1;
                    }
                }
            }
            k2++;
        }
        k1++;
    }

    if (update) {
        if (dbg_lvl()) {
            sol.print_solution();
        }
        updated = true;

        sol.update_insertion_move_inter(i_best, j_best, k_best, kk_best, l);
        if (dbg_lvl()) {
            sol.print_solution();
        }
    } else {
        updated = false;
    }

    if (dbg_lvl() > 0) {
        printf(
            "inter insertion with l = %d, running time = %f and improvement %d "
            "on machines %d and %d on places %d %d\n",
            l, CCutil_zeit() - runningtime, max, k_best, kk_best, i_best,
            j_best);
        print_line();
    }
}
/** update to best improvement */

void LocalSearchData::calculate_W(const Sol& sol) {
    std::size_t i = 0UL;
    for (auto& m : sol.machines) {
        auto tmp = m.job_list.begin();
        W[i][0] = value_Fj(sol.c[(*tmp)->job], (*tmp).get());
        tmp++;
        int j = 1;
        for (; tmp != m.job_list.end(); tmp++) {
            W[i][j] = W[i][j - 1] + value_Fj(sol.c[(*tmp)->job], (*tmp).get());
            j++;
        }
        i++;
    }
}

void LocalSearchData::calculate_g(const Sol& sol) {
    int i = 0;
    for (auto& m : sol.machines) {
        for (int j = 0; j < nb_jobs; j++) {
            g[i][j].clear();
        }  // GPtrArray*  machine1 = sol->part[k1].machine;

        int P{0};
        int j{0};
        for (auto it = m.job_list.begin(); it != m.job_list.end(); it++) {
            VecJobPtr lateness_sort(it, m.job_list.end());

            std::sort(lateness_sort.begin(), lateness_sort.end(),
                      [&sol](const auto& lhs, const auto& rhs) -> bool {
                          return (sol.c[lhs->job] - lhs->due_time >
                                  sol.c[rhs->job] - rhs->due_time);
                      });

            int  tw{0};
            int  w{0};
            int  t1{0};
            bool move{true};
            auto k = 0U;
            auto tmp = lateness_sort[k];

            for (; k < lateness_sort.size() && move;) {
                move = tmp->weight * (sol.c[tmp->job] - P - tmp->due_time) > 0;

                if (move) {
                    tw += tmp->weight * (sol.c[tmp->job] - P - tmp->due_time);
                    w += tmp->weight;
                    k++;
                    tmp = lateness_sort[k];
                }
            }

            int t2 = tmp->due_time - sol.c[tmp->job] + P;

            g[i][j].emplace_back(t1, t2, tw, w);

            for (auto l = k; l < lateness_sort.size();) {
                tw = tw + w * (t2 - t1);
                t1 = t2;
                move = true;
                tmp = lateness_sort[l];

                while (move) {
                    w += tmp->weight;
                    l++;

                    if (l == lateness_sort.size()) {
                        move = false;
                        t2 = std::numeric_limits<int>::max();
                    } else {
                        tmp = lateness_sort[l];
                        t2 = tmp->due_time - sol.c[tmp->job] + P;
                        move = (t1 == t2);
                    }
                }

                g[i][j].emplace_back(t1, t2, tw, w);
            }

            tmp = m.job_list[j];
            P += tmp->processing_time;
            j++;
        }
        i++;
    }
}

void LocalSearchData::create_processing_insertion_operator(const Sol& sol,
                                                           int        l) {
    for (auto i = 0UL; auto& m : sol.machines) {
        auto  C = 0;
        auto& machine = m.job_list;

        for (auto j = 0UL; j < l; j++) {
            C += machine[j]->processing_time;
        }

        for (size_t j = 0; j < machine.size() - l; j++) {
            processing_list[0][i][j] = {j, C};
            C = C - machine[j]->processing_time +
                machine[j + l]->processing_time;
        }

        std::sort(processing_list[0][i].begin(),
                  processing_list[0][i].begin() + machine.size() - l,
                  [](const auto& lhs, const auto& rhs) -> bool {
                      return lhs.p > rhs.p;
                  });
        i++;
    }
}

void LocalSearchData::create_processing_list_swap(const Sol& sol,
                                                  int        l1,
                                                  int        l2) {
    int val = 0;

    for (auto i = 0UL; auto& m : sol.machines) {
        auto&                vec_m = m.job_list;
        std::shared_ptr<Job> j1{}, j2{};
        const auto           n = vec_m.size();
        auto                 C = 0;

        for (std::size_t j = l1; j < l1 + l2; ++j) {
            C += vec_m[j]->processing_time;
        }

        for (std::size_t j = l1; j < n - l2; ++j) {
            processing_list[1][i][j - l1] = {j, C};
            C = C - vec_m[j]->processing_time + vec_m[j + l2]->processing_time;
        }

        processing_list[1][i][n - l2 - l1] = {n - l2, C};
        std::sort(
            processing_list[1][i].begin(),
            processing_list[1][i].begin() + n - l1 - l2 + 1,
            [](const auto& lhs, const auto& rhs) { return lhs.p < rhs.p; });

        C = 0;

        for (int j = 0; j < l1; ++j) {
            C += vec_m[j]->processing_time;
        }

        for (std::size_t j = 0; j < n - l1 - l2; ++j) {
            processing_list[0][i][j] = {j, C};
            C = C - vec_m[j]->processing_time + vec_m[j + l1]->processing_time;
        }

        processing_list[0][i][n - l1 - l2] = {n - l1 - l2, C};
        std::sort(
            processing_list[0][i].begin(),
            processing_list[0][i].begin() + n - l1 - l2 + 1,
            [](const auto& lhs, const auto& rhs) { return lhs.p > rhs.p; });
        i++;
    }
}

void LocalSearchData::create_processing_insertion_operator_inter(const Sol& sol,
                                                                 int        l) {
    int val = 0;

    for (auto i = 0UL; auto& m : sol.machines) {
        const auto& nb_jobs = m.job_list.size();
        const auto& vec_m = m.job_list;
        auto        C = 0;

        for (auto j = 0UL; j < l; ++j) {
            auto j1 = vec_m[j];
            C += j1->processing_time;
        }

        for (auto j = 0UL; j < nb_jobs - l; ++j) {
            auto j1 = vec_m[j];
            auto j2 = vec_m[j + l];
            processing_list[0][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[0][i][nb_jobs - l] = {nb_jobs - l, C};
        std::sort(
            processing_list[0][i].begin(),
            processing_list[0][i].begin() + nb_jobs - l + 1,
            [](const auto& lhs, const auto& rhs) { return lhs.p < rhs.p; });
        i++;
    }
}

void LocalSearchData::create_processing_list_swap_inter(const Sol& sol,
                                                        int        l1,
                                                        int        l2) {
    int val = 0;

    for (int i = 0; i < nb_machines; ++i) {
        auto&      machine = sol.machines[i].job_list;
        const auto nb_jobs = machine.size();
        int        C = 0;

        if (nb_jobs - l1 < 0 || nb_jobs - l2 < 0) {
            continue;
        }

        for (int j = 0; j < l1; ++j) {
            C += machine[j]->processing_time;
        }

        for (auto j = 0UL; j < nb_jobs - l1; ++j) {
            auto j1 = machine[j];
            auto j2 = machine[j + l1];
            processing_list[0][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[0][i][nb_jobs - l1].pos = nb_jobs - l1;
        processing_list[0][i][nb_jobs - l1].p = C;
        std::sort(
            processing_list[0][i].begin(),
            processing_list[0][i].begin() + nb_jobs - l1 + 1,
            [](const auto& lhs, const auto& rhs) { return lhs.p < rhs.p; });
        C = 0;

        for (int j = 0; j < l2; ++j) {
            C += machine[j]->processing_time;
        }

        for (auto j = 0UL; j < nb_jobs - l2; ++j) {
            auto j1 = machine[j];
            auto j2 = machine[j + l2];
            processing_list[1][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[1][i][nb_jobs - l2] = {nb_jobs - l2, C};
        std::sort(
            processing_list[1][i].begin(),
            processing_list[1][i].begin() + nb_jobs - l2 + 1,
            [](const auto& lhs, const auto& rhs) { return lhs.p < rhs.p; });
    }
}
