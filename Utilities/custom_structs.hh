#ifndef CUSTOM_STRUCTS_HH
#define CUSTOM_STRUCTS_HH

struct double_int {
    double value;
    int rank;
    double_int() : value(1e200), rank(-1) {};
    double_int(double v, int r): value(v), rank(r) {};
};

#endif
