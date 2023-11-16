#include "vmath.h"
#include <math.h>

double clamp(double value, double lower, double upper) {
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

double lerp(double start, double end, double amount) {
    return (1 - amount) * start + amount * end;
}

static const vector3_t vector3_zero = { 0 };

vector3_t vector3_add(vector3_t const *u, vector3_t const *v) {
    return (vector3_t){ .a[0] = u->a[0] + v->a[0], .a[1] = u->a[1] + v->a[1], .a[2] = u->a[2] + v->a[2] };
}

vector3_t vector3_subtract(vector3_t const *u, vector3_t const *v) {
    return (vector3_t){ .a[0] = u->a[0] - v->a[0], .a[1] = u->a[1] - v->a[1], .a[2] = u->a[2] - v->a[2] };
}

vector3_t vector3_scale(vector3_t const *v, double s) {
    return (vector3_t){ .a[0] = v->a[0] * s, .a[1] = v->a[1] * s, .a[2] = v->a[2] * s };
}

vector3_t vector3_negate(vector3_t const *v) {
    return (vector3_t){ .a[0] = (-1) * v->a[0], .a[1] = (-1) * v->a[1], .a[2] = (-1) * v->a[2] };
}

double vector3_dot(vector3_t const *u, vector3_t const *v) {
    return u->a[0] * v->a[0] + u->a[1] * v->a[1] + u->a[2] * v->a[2];
}

double vector3_length2(vector3_t const *v) {
    return vector3_dot(v, v);
}

double vector3_length(vector3_t const *v) {
    return sqrt(vector3_dot(v, v));
}

vector3_t vector3_normalize(vector3_t const *v) {
    if (v->a[0] == 0 && v->a[1] == 0 && v->a[2] == 0) {
        return vector3_zero;
    }
    double length = vector3_length(v);
    return (vector3_t){ .a[0] = v->a[0] / length, .a[1] = v->a[1] / length, .a[2] = v->a[2] / length };
}

vector3_t vector3_project(vector3_t const *v, vector3_t const *direction) {
    double factor = vector3_dot(v, direction) / vector3_dot(direction, direction);
    return (vector3_t){ .a[0] = direction->a[0] * factor, .a[1] = direction->a[1] * factor, .a[2] = direction->a[2] * factor };
}

vector3_t vector3_reflect(vector3_t const *v, vector3_t const *direction) {
    double factor = 2 * vector3_dot(v, direction) / vector3_dot(direction, direction);
    return (vector3_t
    ){ .a[0] = direction->a[0] * factor - v->a[0], .a[1] = direction->a[1] * factor - v->a[1], .a[2] = direction->a[2] * factor - v->a[2] };
}

vector3_t vector3_cross(vector3_t const *u, vector3_t const *v) {
    return (vector3_t){
        .a[0] = u->a[1] * v->a[2] - u->a[2] * v->a[1],
        .a[1] = u->a[2] * v->a[0] - u->a[0] * v->a[2],
        .a[2] = u->a[0] * v->a[1] - u->a[1] * v->a[0],
    };
}

vector3_t vector3_transform(vector3_t const *v, matrix_t const *transform) {
    return (vector3_t
    ){ .a[0] = transform->a[0][0] * v->a[0] + transform->a[0][1] * v->a[1] + transform->a[0][2] * v->a[2] + transform->a[0][3],
       .a[1] = transform->a[1][0] * v->a[0] + transform->a[1][1] * v->a[1] + transform->a[1][2] * v->a[2] + transform->a[1][3],
       .a[2] = transform->a[2][0] * v->a[0] + transform->a[2][1] * v->a[1] + transform->a[2][2] * v->a[2] + transform->a[2][3] };
}

vector4_t vector4_scale(vector4_t const *v, double s) {
    return (vector4_t){
        .a[0] = v->a[0] * s,
        .a[1] = v->a[1] * s,
        .a[2] = v->a[2] * s,
        .a[3] = v->a[3] * s,
    };
}

vector4_t vector4_lerp(vector4_t const *from, vector4_t const *to, double amount) {
    return (vector4_t){
        .a[0] = (1 - amount) * from->a[0] + amount * to->a[0],
        .a[1] = (1 - amount) * from->a[1] + amount * to->a[1],
        .a[2] = (1 - amount) * from->a[2] + amount * to->a[2],
        .a[3] = (1 - amount) * from->a[3] + amount * to->a[3],
    };
}

vector4_t vector4_clamp(vector4_t const *v, double lower, double upper) {
    return (vector4_t){
        .a[0] = clamp(v->a[0], lower, upper),
        .a[1] = clamp(v->a[1], lower, upper),
        .a[2] = clamp(v->a[2], lower, upper),
        .a[3] = clamp(v->a[3], lower, upper),
    };
}

matrix_t matrix_multiply(matrix_t const *p, matrix_t const *q) {
    matrix_t product = { 0 };
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            product.a[i][j] = 0;
            for (size_t k = 0; k < 4; ++k) {
                product.a[i][j] += p->a[i][k] * q->a[k][j];
            }
        }
    }
    return product;
}

// matrix_t matrix_look_at(vector3_t const *from, vector3_t const *to, vector3_t const *up) {}

// "counter-clockwise" in right-handed coordinate systems
matrix_t matrix_x_rotation_from_angle(double theta) {
    return (matrix_t){
        .a[0][0] = 1,
        .a[0][1] = 0,
        .a[0][2] = 0,
        .a[0][3] = 0,
        .a[1][0] = 0,
        .a[1][1] = cos(theta),
        .a[1][2] = -sin(theta),
        .a[1][3] = 0,
        .a[2][0] = 0,
        .a[2][1] = sin(theta),
        .a[2][2] = cos(theta),
        .a[2][3] = 0,
        .a[3][0] = 0,
        .a[3][1] = 0,
        .a[3][2] = 0,
        .a[3][3] = 1,
    };
}

matrix_t matrix_y_rotation_from_angle(double theta) {
    return (matrix_t){
        .a[0][0] = cos(theta),
        .a[0][1] = 0,
        .a[0][2] = sin(theta),
        .a[0][3] = 0,
        .a[1][0] = 0,
        .a[1][1] = 1,
        .a[1][2] = 0,
        .a[1][3] = 0,
        .a[2][0] = -sin(theta),
        .a[2][1] = 0,
        .a[2][2] = cos(theta),
        .a[2][3] = 0,
        .a[3][0] = 0,
        .a[3][1] = 0,
        .a[3][2] = 0,
        .a[3][3] = 1,
    };
}

matrix_t matrix_z_rotation_from_angle(double theta) {
    return (matrix_t){
        .a[0][0] = cos(theta),
        .a[0][1] = -sin(theta),
        .a[0][2] = 0,
        .a[0][3] = 0,
        .a[1][0] = sin(theta),
        .a[1][1] = cos(theta),
        .a[1][2] = 0,
        .a[1][3] = 0,
        .a[2][0] = 0,
        .a[2][1] = 0,
        .a[2][2] = 1,
        .a[2][3] = 0,
        .a[3][0] = 0,
        .a[3][1] = 0,
        .a[3][2] = 0,
        .a[3][3] = 1,
    };
}
