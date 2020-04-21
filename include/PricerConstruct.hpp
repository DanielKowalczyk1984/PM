#include <glib.h>
#include <interval.h>
#include <solution.h>
#include <boost/dynamic_bitset.hpp>
#include "NodeBddSpec.hpp"

class conflict_state {
   public:
    boost::dynamic_bitset<> add;
    boost::dynamic_bitset<> remove;

    conflict_state(){};

    ~conflict_state(){};
};

class PricerConstruct : public DdSpec<PricerConstruct, int, 2> {
    GPtrArray* pair_list;
    int        nb_layers;

   public:
    explicit PricerConstruct(GPtrArray* _pair_list) : pair_list(_pair_list) {
        nb_layers = pair_list->len;
    };

    int getRoot(int& state) {
        state = 0;
        return nb_layers;
    };

    int getChild(int& state, int level, int value) const {
        int layer = nb_layers - level;
        int _j;
        assert(0 <= layer && layer <= nb_layers - 1);
        job_interval_pair* tmp_pair =
            (job_interval_pair*)g_ptr_array_index(pair_list, layer);
        interval* tmp_interval = tmp_pair->I;
        Job*      tmp_j = (Job*)tmp_pair->j;

        if (value) {
            state = state + tmp_j->processing_time;
        }

        _j = min_job(layer, state, value);

        if (!(_j < nb_layers)) {
            if (state <= tmp_interval->b) {
                return -1;
            }
            return 0;
        }
        assert(_j < nb_layers);
        return nb_layers - _j;
    }

    ~PricerConstruct(){};

   private:
    int min_job(int j, int state, int value) const {
        int  val = nb_layers;
        Job* tmp = ((job_interval_pair*)g_ptr_array_index(pair_list, j))->j;

        if (value) {
            for (int i = j + 1; i < nb_layers; ++i) {
                job_interval_pair* tmp_pair =
                    (job_interval_pair*)g_ptr_array_index(pair_list, i);
                interval* tmp_interval = tmp_pair->I;
                Job*      tmp_j = tmp_pair->j;

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b) {
                    if (tmp == tmp_j ) {
                        continue;
                    }
                    val = i;
                    break;
                }
            }
        } else {
            for (int i = j + 1; i < nb_layers; ++i) {
                job_interval_pair* tmp_pair =
                    (job_interval_pair*)g_ptr_array_index(pair_list, i);
                interval* tmp_interval = tmp_pair->I;
                Job*      tmp_j = tmp_pair->j;

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b) {
                    val = i;
                    break;
                }
            }
        }

        return val;
    }

    int diff_obj(Job* i, Job* j, int C) const {
        return value_Fj(C, i) + value_Fj(C + j->processing_time, j) -
               (value_Fj(C - i->processing_time + j->processing_time, j) +
                value_Fj(C + j->processing_time, i));
    }
};

class PricerConstructTI : public DdSpec<PricerConstructTI, int, 2> {
    GPtrArray* pair_list;
    int*       take_job;
    int        Hmax;
    int        nb_layers;

   public:
    explicit PricerConstructTI(GPtrArray* _pair_list, int* _take_job, int _Hmax)
        : pair_list(_pair_list), take_job(_take_job), Hmax(_Hmax) {
        nb_layers = pair_list->len;
    };

    int getRoot(int& state) {
        state = 0;
        return nb_layers;
    };

    int getChild(int& state, int level, int value) const {
        int layer = nb_layers - level;
        int _j;
        assert(0 <= layer && layer <= nb_layers - 1);
        job_interval_pair* tmp_pair =
            (job_interval_pair*)g_ptr_array_index(pair_list, layer);
        interval* tmp_interval = tmp_pair->I;
        Job*      tmp_j = (Job*)tmp_pair->j;

        if (level - 1 == 0 && value) {
            return (state + tmp_j->processing_time <= tmp_interval->b)? -1 :
            0;
        } else if (level - 1 == 0) {
            return ( state <= tmp_interval->b) ? -1 : 0;
        }

        if (value) {
            state = state + tmp_j->processing_time;
        }

        _j = min_job(layer, state, value);

        if (!(_j < nb_layers)) {
            if (state <= tmp_interval->b) {
                return -1;
            }
            return 0;
        }
        assert(_j < nb_layers);
        return nb_layers - _j;
    }

    ~PricerConstructTI(){};

   private:
    int min_job(int j, int state, int value) const {
        int                val = nb_layers;
        job_interval_pair* tmp_pair;
        interval*          tmp_interval;
        Job*               tmp_j;
        Job* tmp = ((job_interval_pair*)g_ptr_array_index(pair_list, j))->j;

        if (value) {
            for (int i = j + 1; i < nb_layers; ++i) {
                tmp_pair = (job_interval_pair*)g_ptr_array_index(pair_list, i);
                tmp_interval = tmp_pair->I;
                tmp_j = tmp_pair->j;

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b &&
                    take_job[tmp_j->job * (Hmax + 1) + state]) {
                    if (tmp == tmp_j ||
                        (tmp->job > tmp_j->job &&
                         value_diff_Fij(state, tmp_j, tmp) <= 0)) {
                        continue;
                    }
                    val = i;
                    break;
                }
            }
        } else {
            for (int i = j + 1; i < nb_layers; ++i) {
                tmp_pair = (job_interval_pair*)g_ptr_array_index(pair_list, i);
                tmp_interval = tmp_pair->I;
                tmp_j = tmp_pair->j;

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b &&
                    take_job[tmp_j->job * (Hmax + 1) + state]) {
                    val = i;
                    break;
                }
            }
        }

        return val;
    }

    int diff_obj(Job* i, Job* j, int C) const {
        return value_Fj(C, i) + value_Fj(C + j->processing_time, j) -
               (value_Fj(C - i->processing_time + j->processing_time, j) +
                value_Fj(C + j->processing_time, i));
    }
};

class ConflictConstraints
    : public DdSpec<ConflictConstraints, conflict_state, 2> {
    int                                  nb_jobs;
    std::vector<boost::dynamic_bitset<>> differsets;
    std::vector<boost::dynamic_bitset<>> samesets;

    bool takeable(int job, conflict_state& state) {
        if (state.remove[job]) {
            return false;
        }

        return true;
    }

    bool leaveable(int job, conflict_state& state) {
        if (state.add[job]) {
            return false;
        }

        return true;
    }

   public:
    ConflictConstraints(int _nb_jobs, int* elist_same, int edge_count_same,
                        int* elist_differ, int edge_count_differ)
        : nb_jobs(_nb_jobs) {
        differsets.resize(_nb_jobs);
        samesets.resize(_nb_jobs);

        for (int i = 0; i < _nb_jobs; i++) {
            differsets[i].resize(_nb_jobs);
            samesets[i].resize(_nb_jobs);
        }

        for (int i = 0; i < edge_count_same; ++i) {
            samesets[elist_same[2 * i]][elist_same[2 * i + 1]] = 1;
        }

        for (int i = 0; i < edge_count_differ; ++i) {
            differsets[elist_differ[2 * i]][elist_differ[2 * i + 1]] = 1;
        }
    };

    ~ConflictConstraints() {}

    int getRoot(conflict_state& state) const {
        state.add.resize(nb_jobs);
        state.remove.resize(nb_jobs);
        return nb_jobs;
    }

    int getChild(conflict_state& state, int level, int take) {
        int job = nb_jobs - level;
        int _j;
        assert(0 <= job && job <= nb_jobs - 1);

        if (samesets[job].intersects(differsets[job])) {
            return 0;
        }

        if (level - 1 == 0 && take) {
            return (!state.remove[job]) ? -1 : 0;
        } else if (level - 1 == 0) {
            return (!state.add[job]) ? -1 : 0;
        }

        if (take) {
            if (!takeable(job, state)) {
                return 0;
            }

            state.add |= samesets[job];
            state.remove |= differsets[job];
        } else {
            if (!leaveable(job, state)) {
                return 0;
            }

            state.remove |= samesets[job];
        }

        _j = min_job(job, state);

        if (_j == nb_jobs && take) {
            return (!state.remove[job]) ? -1 : 0;
        } else if (_j == nb_jobs) {
            return (!state.add[job]) ? -1 : 0;
        }

        assert(_j < nb_jobs);
        return nb_jobs - _j;
    }

    bool equalTo(conflict_state const& state1,
                 conflict_state const& state2) const {
        if (state2.add != state1.add) {
            return false;
        }

        if (state2.remove != state1.remove) {
            return false;
        }

        return true;
    }

    size_t hashCode(conflict_state const& state) const {
        size_t val = 0;
        size_t it = state.add.find_first();

        while (it != boost::dynamic_bitset<>::npos) {
            val += 1213657 * static_cast<size_t>(it);
            it = state.add.find_next(it);
        }

        it = state.remove.find_first();

        while (it != boost::dynamic_bitset<>::npos) {
            val += 487239 * static_cast<size_t>(it);
            it = state.remove.find_next(it);
        }

        return val;
    }

    int min_job(int j, conflict_state& state) const {
        int i, val = nb_jobs;

        for (i = j + 1; i < nb_jobs; ++i) {
            if (!state.remove[i]) {
                val = i;
                break;
            }
        }

        return val;
    }
};

class scheduling : public DdSpec<scheduling, int, 2> {
    Job*       job;
    GPtrArray* list_layers;
    int        nb_layers;
    int        order;

   public:
    scheduling(Job* _job, GPtrArray* _list_layers, int _order)
        : job(_job), list_layers(_list_layers), order(_order) {
        nb_layers = list_layers->len;
    };

    ~scheduling() {}

    int getRoot(int& state) const {
        state = 0;
        return nb_layers;
    }

    int getChild(int& state, int level, int take) {
        int                j = nb_layers - level;
        job_interval_pair* tmp =
            (job_interval_pair*)g_ptr_array_index(list_layers, j);
        Job* tmp_j = tmp->j;
        assert(0 <= j && j <= nb_layers - 1);

        if (level - 1 == 0 && take) {
            if (tmp_j != job) {
                return -1;
            } else {
                if (state == 0) {
                    return -1;
                }
                return 0;
            }
        } else if (level - 1 == 0) {
            return -1;
        }

        if (take) {
            if (tmp_j != job) {
                // state++;
                j++;
                return nb_layers - j;
            } else {
                if (state == 0) {
                    state++;
                    j++;
                    return nb_layers - j;
                } else if (state >= order) {
                    state = 1;
                    j++;
                    return nb_layers - j;
                } else {
                    state++;
                    return 0;
                }
            }
        } else {
            j++;
            return nb_layers - j;
        }
    }
};
