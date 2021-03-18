#include "Job.h"
#include "util.h"

Job::Job(int p, int w, int d) : processing_time(p), due_time(d), weight(w) {}

extern inline int value_Fj(int C, Job* j);

int value_diff_Fij(int C, Job* i, Job* j) {
    auto val = value_Fj(C + i->processing_time - j->processing_time, i);
    val += value_Fj(C + i->processing_time, j);
    val -= value_Fj(C, j);
    val -= value_Fj(C + i->processing_time, i);
    return val;
}

int bool_diff_Fij(int weight, Job* _prev, Job* tmp_j) {
    return (_prev == nullptr) ? 1
                              : (value_diff_Fij(weight + tmp_j->processing_time,
                                                _prev, tmp_j) >= 0);
}

int arctime_diff_Fij(int weight, Job* i, Job* j) {
    auto val = value_Fj(weight, i);
    val += value_Fj(weight + j->processing_time, j);
    val -= value_Fj(weight - i->processing_time + j->processing_time, j);
    val -= value_Fj(weight + j->processing_time, i);
    return val;
}
