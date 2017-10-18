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
    int nlayers;

public:
    PricerConstruct(GPtrArray *_pair_list)
    : pair_list(_pair_list){
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
             return (state + tmp_j->processingime <= tmp_interval->b)? -1 : 0;
         } else if (level - 1 == 0) {
             return ( state <= tmp_interval->b) ? -1 : 0;
         }

         if (value) {
             state = state + tmp_j->processingime;
             _j = min_job(layer, state, value);

             if(!(_j < nlayers)) {
                if( state <= tmp_interval->b) {
                    if(value_Fj(state, tmp_j) - value_Fj(state + 1, tmp_j) > 0) {
                        return 0;
                    }
                    return -1;
                }
                return 0;
             }
         } else {
             _j = min_job(layer, state,value);

             if(!(_j < nlayers)) {
                if( state <= tmp_interval->b) {
                    return -1;
                }
                return 0;
             }
         }

         assert(_j < nlayers);
         return nlayers - _j;
     }

    ~PricerConstruct(){};

    private:
     int min_job(int j, int state, int value) const {
         int  val = nlayers;
         job_interval_pair *tmp_pair;
         interval *tmp_interval;
         Job *tmp_j;
         Job * tmp = ((job_interval_pair*) g_ptr_array_index(pair_list,j))->j;
         int key = ((job_interval_pair*) g_ptr_array_index(pair_list,j))->I->key;

         if(value) {
             for (int i = j + 1; i < nlayers; ++i) {
                tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list, i);
                tmp_interval = tmp_pair->I;
                tmp_j = tmp_pair->j;

                 if (state + tmp_j->processingime > tmp_interval->a  && state + tmp_j->processingime <= tmp_interval->b ) {
                     if(diff_obj(tmp, tmp_j, state) >= 0 && key != tmp_interval->key) {
                         continue;
                     }
                     val = i;
                     break;
                 }
             }
         } else {
            for (int i = j + 1; i < nlayers; ++i) {
               tmp_pair = (job_interval_pair *) g_ptr_array_index(pair_list, i);
               tmp_interval = tmp_pair->I;
               tmp_j = tmp_pair->j;

                if (state + tmp_j->processingime > tmp_interval->a  && state + tmp_j->processingime <= tmp_interval->b) {
                    val = i;
                    break;
                }
            }
         }

         return val;
     }

     int diff_obj(Job *i, Job *j, int C) const {
        return value_Fj(C, i) + value_Fj(C + j->processingime, j) - (value_Fj(C - i->processingime + j->processingime, j) + value_Fj(C + j->processingime, i));

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

class scheduling: public tdzdd::DdSpec<scheduling, int, 2> {
    Job * job;
    GPtrArray *list_layers;
    int nlayers;
    int order;


   public:
    scheduling(Job *_job, GPtrArray *_list_layers, int _order):job(_job) ,list_layers(_list_layers),order(_order) {
        nlayers = list_layers->len;
    };

    ~scheduling() {}

    int getRoot(int &state) const {
        state = 0;
        return nlayers;
    }

    int getChild(int &state, int level, int take) {
        int j = nlayers - level;
        job_interval_pair *tmp = (job_interval_pair *) g_ptr_array_index(list_layers, j);
        Job *tmp_j = tmp->j;
        assert(0 <= j && j <= nlayers - 1);


        if (level - 1 == 0 && take) {
            if(tmp_j != job) {
                return -1;
            } else {
                if(state == 0) {
                    return -1;
                }
                return 0;
            }
        } else if (level - 1 == 0) {
            return -1;
        }

        if (take) {
            if(tmp_j != job) {
                j++;
                return nlayers - j;
            } else {
                if(state == 0) {
                    state++;
                    j++;
                    return  nlayers - j;
                }else if (state < order) {
                    state++;
                    j++;
                    return nlayers - j;

                } else {
                    return 0;
                }
            }
        } else {
            j++;
            return nlayers - j;
        }
    }
};
