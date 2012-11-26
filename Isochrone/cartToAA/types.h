#ifndef __TYPES_H__
#define __TYPES_H__

struct body_struct {
    size_t id;
    double m;
    double q[3];
    double p[3];
    double temp[3];
};
typedef struct body_struct body;

#endif //__TYPES_H__