#include <boost/dynamic_bitset.hpp>
#include <tdzdd/DdSpec.hpp>
#include <vector>
#include "solution.h"
#include "interval.h"
#include <glib.h>

class conflict_state {
   public:
    boost::dynamic_bitset<> add;
    boost::dynamic_bitset<> remove;

    conflict_state(){};

    ~conflict_state(){};
};

class PricerConstruct : public tdzdd::DdSpec<PricerConstruct, int, 2> {
    GPtrArray *pair_list;
    int **sum_p;
    int nlayers;

public:
    PricerConstruct(GPtrArray *_pair_list, int **_sum_p)
    : pair_list(_pair_list), sum_p(_sum_p){
        nlayers = pair_list->len;
    };

    int getRoot(int &state){
        state = 0;
        return nlayers;
    };

     int getChild(int &state, int level, int value) const {
         int layer = nlayers - level;
         int _j;
         assert(0 <= layer && layer <= nlayers - 1);
         job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list, layer);
         interval *tmp_interval = tmp_pair->I;
         Job *tmp_j = (Job *) tmp_pair->j;

         if (level - 1 == 0 && value) {
             return (state + tmp_j->processingime> tmp_interval->a &&
                     state + tmp_j->processingime<= tmp_interval->b)
                        ? -1
                        : 0;
         } else if (level - 1 == 0) {
             return (state <= tmp_interval->b) ? -1 : 0;
         }

         if (value) {
             int sum = state + tmp_j->processingime;
             _j = min_job(layer, sum);

             if(_j < nlayers) {
                 tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list,_j);
                 tmp_interval = tmp_pair->I;
                 tmp_j = tmp_pair->j;
                 if(sum + sum_p[tmp_interval->key][tmp_j->job] < tmp_interval->b) {
                     return 0;
                 }
             }

             if (!(_j < nlayers)) {
                tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list,_j - 1);
                tmp_interval = tmp_pair->I;
                 if ((sum <= tmp_interval->b)) {
                     return -1;
                 }

                 return 0;
             }
             state = sum;
         } else {
             _j = min_job(layer, state);

             if(_j < nlayers) {
                tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list,_j);
                tmp_interval = tmp_pair->I;
                tmp_j = tmp_pair->j;
                if(state + sum_p[tmp_interval->key][tmp_j->job] < tmp_interval->b) {
                    return 0;
                }
             }

             if (!(_j < nlayers)) {
                tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list,_j - 1);
                tmp_interval = tmp_pair->I;
                 if (state <= tmp_interval->b) {
                     return -1;
                 }

                 return 0;
             }
         }

         // if (_j == nbjobs && state >= Hmin && state <= Hmax) {
         //     return -1;
         // } else if (_j == nbjobs) {
         //     return 0;
         // }

         assert(_j < nlayers);
         return nlayers - _j;
     }

    ~PricerConstruct(){};

    private:
     int min_job(int j, int state) const {
         int  val = j + 1;
         job_interval_pair *tmp_pair;
         interval *tmp_interval;
         tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list, j);
         Job *cur_job = tmp_pair->j;
         Job *job;


         for (int i = j + 1; i < nlayers; ++i) {
            tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list, i);
            tmp_interval = tmp_pair->I;
            job = tmp_pair->j;

             if (state > tmp_interval->a - job->processingime &&
                 state <= tmp_interval->b - job->processingime && cur_job->job != job->job) {
                 val = i;
                 break;
             }
         }

         return val;
     }


};

class PricerSpec : public tdzdd::DdSpec<PricerSpec, int, 2> {
    Job *jobarray;
    int  nbjobs;
    int  Hmin;
    int  Hmax;

   public:
    int *sum_p;
    int *min_p;
    PricerSpec(Job *_jobarray, int _nbjobs, int Hmin, int Hmax)
        : jobarray(_jobarray), nbjobs(_nbjobs), Hmin(Hmin), Hmax(Hmax) {
        nbjobs = _nbjobs;
        sum_p = new int[nbjobs];
        min_p = new int[nbjobs];
        int end = nbjobs - 1;
        sum_p[end] = jobarray[end].processingime;
        min_p[end] = jobarray[end].processingime;

        for (int i = end - 1; i >= 0; i--) {
            sum_p[i] = sum_p[i + 1] + jobarray[i].processingime;

            if (jobarray[i].processingime < min_p[i + 1]) {
                min_p[i] = jobarray[i].processingime;
            } else {
                min_p[i] = min_p[i + 1];
            }
        }
    }

    ~PricerSpec() {
        delete[] min_p;
        delete[] sum_p;
    }

    int getRoot(int &state) const {
        state = 0;
        return nbjobs;
    }

    int getChild(int &state, int level, int value) const {
        int job = nbjobs - level;
        int _j;
        assert(0 <= job && job <= nbjobs - 1);
        int temp_p = jobarray[job].processingime;

        if (level - 1 == 0 && value) {
            return (state + jobarray[job].processingime >= Hmin &&
                    state + jobarray[job].processingime <= Hmax)
                       ? -1
                       : 0;
        } else if (level - 1 == 0) {
            return (state >= Hmin && state <= Hmax) ? -1 : 0;
        }

        if (value) {
            int sum = state + temp_p;
            _j = min_job(job, sum);

            if (_j < nbjobs) {
                if (state + sum_p[_j] < Hmin) {
                    return 0;
                }

                if ((sum >= Hmin && sum <= Hmax) && (sum + min_p[_j] > Hmax)) {
                    return -1;
                }
            } else {
                if ((sum >= Hmin && sum <= Hmax)) {
                    return -1;
                }

                return 0;
            }

            state = sum;
        } else {
            _j = min_job(job, state);

            if (_j < nbjobs) {
                if (state + sum_p[_j] < Hmin) {
                    return 0;
                }

                if ((state >= Hmin && state <= Hmax) &&
                    (state + min_p[_j] > Hmax)) {
                    return -1;
                }
            } else {
                if (state >= Hmin && state <= Hmax) {
                    return -1;
                }

                return 0;
            }
        }

        // if (_j == nbjobs && state >= Hmin && state <= Hmax) {
        //     return -1;
        // } else if (_j == nbjobs) {
        //     return 0;
        // }

        assert(_j < nbjobs);
        return nbjobs - _j;
    }

   private:
    int min_job(int j, int state) const {
        int  val = j + 1;

        // for (i = j + 1; i < nbjobs; ++i) {
        //     if (state >= jobarray[i].releasetime &&
        //         state <= jobarray[i].duetime - jobarray[i].processingime) {
        //         val = i;
        //         break;
        //     }
        // }

        return val;
    }
};

class ConflictConstraints
    : public tdzdd::DdSpec<ConflictConstraints, conflict_state, 2> {
    int                                  nbjobs;
    std::vector<boost::dynamic_bitset<>> differsets;
    std::vector<boost::dynamic_bitset<>> samesets;

    bool takeable(int job, conflict_state &state) {
        if (state.remove[job]) {
            return false;
        }

        return true;
    }

    bool leaveable(int job, conflict_state &state) {
        if (state.add[job]) {
            return false;
        }

        return true;
    }

   public:
    ConflictConstraints(int  _nbjobs,
                        int *elist_same,
                        int  ecount_same,
                        int *elist_differ,
                        int  ecount_differ)
        : nbjobs(_nbjobs) {
        differsets.resize(_nbjobs);
        samesets.resize(_nbjobs);

        for (int i = 0; i < _nbjobs; i++) {
            differsets[i].resize(_nbjobs);
            samesets[i].resize(_nbjobs);
        }

        for (int i = 0; i < ecount_same; ++i) {
            samesets[elist_same[2 * i]][elist_same[2 * i + 1]] = 1;
        }

        for (int i = 0; i < ecount_differ; ++i) {
            differsets[elist_differ[2 * i]][elist_differ[2 * i + 1]] = 1;
        }
    };

    ~ConflictConstraints() {}

    int getRoot(conflict_state &state) const {
        state.add.resize(nbjobs);
        state.remove.resize(nbjobs);
        return nbjobs;
    }

    int getChild(conflict_state &state, int level, int take) {
        int job = nbjobs - level;
        int _j;
        assert(0 <= job && job <= nbjobs - 1);

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

        if (_j == nbjobs && take) {
            return (!state.remove[job]) ? -1 : 0;
        } else if (_j == nbjobs) {
            return (!state.add[job]) ? -1 : 0;
        }

        assert(_j < nbjobs);
        return nbjobs - _j;
    }

    bool equalTo(conflict_state const &state1,
                 conflict_state const &state2) const {
        if (state2.add != state1.add) {
            return false;
        }

        if (state2.remove != state1.remove) {
            return false;
        }

        return true;
    }

    size_t hashCode(conflict_state const &state) const {
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

    int min_job(int j, conflict_state &state) const {
        int i, val = nbjobs;

        for (i = j + 1; i < nbjobs; ++i) {
            if (!state.remove[i]) {
                val = i;
                break;
            }
        }

        return val;
    }
};
