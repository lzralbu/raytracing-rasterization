#ifndef __VMATH_H__
#define __VMATH_H__

typedef struct vector2_t {
    double a[2];
} vector2_t;

typedef struct vector3_t {
    double a[3];
} vector3_t;

typedef struct vector4_t {
    double a[4];
} vector4_t;

typedef vector4_t color_t;
typedef vector4_t quaternion_t;

typedef struct matrix_t {
    double a[4][4];
} matrix_t;

double clamp(double value, double lower, double upper);
double lerp(double start, double end, double amount);

void vector3_add(vector3_t *result, vector3_t const *u, vector3_t const *v);
void vector3_subtract(vector3_t *result, vector3_t const *u, vector3_t const *v);

void vector3_scale(vector3_t *result, vector3_t const *v, double s);
void vector3_negate(vector3_t *result, vector3_t const *v);

double vector3_dot(vector3_t const *u, vector3_t const *v);
double vector3_length2(vector3_t const *v);
double vector3_length(vector3_t const *v);
void vector3_normalize(vector3_t *result, vector3_t const *v);

void vector3_project(vector3_t *result, vector3_t const *v, vector3_t const *direction);
void vector3_reflect(vector3_t *result, vector3_t const *v, vector3_t const *direction);

void vector3_cross(vector3_t *result, vector3_t const *u, vector3_t const *v);

void vector3_transform(vector3_t *result, vector3_t const *v, matrix_t const *transform);

void vector4_scale(vector4_t *result, vector4_t const *v, double s);
void vector4_lerp(vector4_t *result, vector4_t const *from, vector4_t const *to, double amount);
void vector4_clamp(vector4_t *result, vector4_t const *v, double lower, double upper);

void matrix_multiply(matrix_t *result, matrix_t const *p, matrix_t const *q);

// void matrix_look_at(matrix_t *result, vector3_t const *from, vector3_t const *to, vector3_t const *up);

void matrix_x_rotation_from_angle(matrix_t *result, double theta);
void matrix_y_rotation_from_angle(matrix_t *result, double theta);
void matrix_z_rotation_from_angle(matrix_t *result, double theta);

#endif // __VMATH_H__
