#include <WeightZDD.hpp>
#include <FarkasZDD.hpp>
#include <DurationBDD.hpp>




struct WeightZDDdouble : WeightZDD<WeightZDDdouble, double> {
    WeightZDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : WeightZDD<WeightZDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};

struct FarkasZDDdouble : FarkasZDD<FarkasZDDdouble, double> {
    FarkasZDDdouble(
        double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : FarkasZDD<FarkasZDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};

struct DurationBDDdouble : DurationBDD<DurationBDDdouble, double> {
    DurationBDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : DurationBDD<DurationBDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};

