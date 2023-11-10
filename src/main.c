#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "raylib.h"

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

typedef struct vector2_t {
    double x;
    double y;
} vector2_t;

typedef struct matrix3_t {
    double a11;
    double a12;
    double a13;
    double a21;
    double a22;
    double a23;
    double a31;
    double a32;
    double a33;
} matrix3_t;

matrix3_t matrix3_multiply(matrix3_t p, matrix3_t q) {
    return (matrix3_t){
        .a11 = p.a11 * q.a11 + p.a12 * q.a21 + p.a13 * q.a31,
        .a12 = p.a11 * q.a12 + p.a12 * q.a22 + p.a13 * q.a32,
        .a13 = p.a11 * q.a13 + p.a12 * q.a23 + p.a13 * q.a33,
        .a21 = p.a21 * q.a11 + p.a22 * q.a21 + p.a23 * q.a31,
        .a22 = p.a21 * q.a12 + p.a22 * q.a22 + p.a23 * q.a32,
        .a23 = p.a21 * q.a13 + p.a22 * q.a23 + p.a23 * q.a33,
        .a31 = p.a31 * q.a11 + p.a32 * q.a21 + p.a33 * q.a31,
        .a32 = p.a31 * q.a12 + p.a32 * q.a22 + p.a33 * q.a32,
        .a33 = p.a31 * q.a13 + p.a32 * q.a23 + p.a33 * q.a33,
    };
}

// "counter-clockwise" in right-handed coordinate systems
matrix3_t matrix3_x_rotation_from_angle(double theta) {
    return (matrix3_t){
        .a11 = 1,
        .a12 = 0,
        .a13 = 0,
        .a21 = 0,
        .a22 = cos(theta),
        .a23 = -sin(theta),
        .a31 = 0,
        .a32 = sin(theta),
        .a33 = cos(theta),
    };
}

matrix3_t matrix3_y_rotation_from_angle(double theta) {
    return (matrix3_t){
        .a11 = cos(theta),
        .a12 = 0,
        .a13 = sin(theta),
        .a21 = 0,
        .a22 = 1,
        .a23 = 0,
        .a31 = -sin(theta),
        .a32 = 0,
        .a33 = cos(theta),
    };
}

matrix3_t matrix3_z_rotation_from_angle(double theta) {
    return (matrix3_t){
        .a11 = cos(theta),
        .a12 = -sin(theta),
        .a13 = 0,
        .a21 = sin(theta),
        .a22 = cos(theta),
        .a23 = 0,
        .a31 = 0,
        .a32 = 0,
        .a33 = 1,
    };
}

typedef struct vector3_t {
    double x;
    double y;
    double z;
} vector3_t;

vector3_t vector3_add(vector3_t v1, vector3_t v2) {
    return (vector3_t){ .x = v1.x + v2.x, .y = v1.y + v2.y, .z = v1.z + v2.z };
}

vector3_t vector3_subtract(vector3_t v1, vector3_t v2) {
    return (vector3_t){ .x = v1.x - v2.x, .y = v1.y - v2.y, .z = v1.z - v2.z };
}

vector3_t vector3_scale(vector3_t v, double s) {
    return (vector3_t){ .x = v.x * s, .y = v.y * s, .z = v.z * s };
}

vector3_t vector3_negate(vector3_t v) {
    return vector3_scale(v, -1.0);
}

double vector3_dot(vector3_t v1, vector3_t v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double vector3_length(vector3_t v) {
    return sqrt(vector3_dot(v, v));
}

vector3_t vector3_normalize(vector3_t v) {
    return vector3_scale(v, 1 / vector3_length(v));
}

vector3_t vector3_reflect(vector3_t v, vector3_t normal) {
    return vector3_subtract(vector3_scale(normal, 2 * vector3_dot(v, normal)), v);
}

vector3_t vector3_transform(vector3_t v, matrix3_t transform) {
    return (vector3_t){
        .x = transform.a11 * v.x + transform.a12 * v.y + transform.a13 * v.z,
        .y = transform.a21 * v.x + transform.a22 * v.y + transform.a23 * v.z,
        .z = transform.a31 * v.x + transform.a32 * v.y + transform.a33 * v.z
    };
}

typedef struct camera_t {
    matrix3_t rotation;
    vector3_t position;
} camera_t;

typedef struct color_t {
    double r;
    double g;
    double b;
} color_t;

typedef struct ray_t {
    vector3_t position;
    vector3_t direction;
} ray_t;

typedef struct sphere_t {
    vector3_t center;
    double radius;
    color_t color;
    double specular;
    double reflective;
} sphere_t;

typedef struct ray_sphere_intersection_t {
    sphere_t const *sphere;
    double t;
} ray_sphere_intersection_t;

enum LIGHT_TYPE {
    LIGHT_TYPE_POINT,
    LIGHT_TYPE_DIRECTIONAL,
    LIGHT_TYPE_AMBIENT
};

typedef struct light_t {
    int type;
    double intensity;
    union {
        vector3_t position;
        vector3_t direction;
    } data;
} light_t;

static const int CANVAS_WIDTH = 900;
static const int CANVAS_HEIGHT = 900;

static const double VIEWPORT_WIDTH = 1;
static const double VIEWPORT_HEIGHT = 1;
static const double VIEWPORT_CAMERA_DISTANCE = 1;

#define BACKGROUND_COLOR \
    (color_t) { 0, 0, 0 }

static sphere_t spheres[] = {
    [0] = { .center = { 0, -1, 3 }, .radius = 1, .color = { 1, 0, 0 }, .specular = 500, .reflective = 0.2f },
    [1] = { .center = { 2, 0, 4 }, .radius = 1, .color = { 0, 0, 1 }, .specular = 500, .reflective = 0.3f },
    [2] = { .center = { -2, 0, 4 }, .radius = 1, .color = { 0, 1, 0 }, .specular = 10, .reflective = 0.4f },
    [3] = { .center = { 0, -5001, 0 }, .radius = 5000, .color = { 1, 1, 0 }, .specular = 1000, .reflective = 0.5f }
};
#define SPHERES_SIZE (sizeof(spheres) / sizeof(sphere_t))

static light_t lights[] = {
    [0] = { .type = LIGHT_TYPE_AMBIENT, .intensity = 0.2 },
    [1] = { .type = LIGHT_TYPE_POINT, .intensity = 0.6, .data.position = { 2, 1, 0 } },
    [2] = { .type = LIGHT_TYPE_DIRECTIONAL, .intensity = 0.2, .data.direction = { 1, 4, 4 } }
};
#define LIGHTS_SIZE (sizeof(lights) / sizeof(light_t))

vector2_t intersect_ray_sphere(ray_t ray, sphere_t sphere) {
    vector3_t temp_vec = vector3_subtract(ray.position, sphere.center);

    double a = vector3_dot(ray.direction, ray.direction);
    double b = 2 * vector3_dot(temp_vec, ray.direction);
    double c = vector3_dot(temp_vec, temp_vec) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return (vector2_t){ .x = INFINITY, .y = INFINITY };
    }

    double t1 = (-b + sqrt(discriminant)) / (2.0 * a);
    double t2 = (-b - sqrt(discriminant)) / (2.0 * a);
    return (vector2_t){ .x = t1, .y = t2 };
}

ray_sphere_intersection_t closest_intersection(ray_t ray, double t_min, double t_max) {
    double closest_t = INFINITY;
    sphere_t const *closest_sphere = 0;
    for (size_t i = 0; i < SPHERES_SIZE; ++i) {
        vector2_t intersection_parameters = intersect_ray_sphere(ray, spheres[i]);
        if (t_min <= intersection_parameters.x && intersection_parameters.x <= t_max && intersection_parameters.x < closest_t) {
            closest_t = intersection_parameters.x;
            closest_sphere = spheres + i;
        }

        if (t_min <= intersection_parameters.y && intersection_parameters.y <= t_max && intersection_parameters.y < closest_t) {
            closest_t = intersection_parameters.y;
            closest_sphere = spheres + i;
        }
    }
    return (ray_sphere_intersection_t){ .sphere = closest_sphere, .t = closest_t };
}

bool is_shadowed(ray_t ray, double t_min, double t_max) {
    for (size_t i = 0; i < SPHERES_SIZE; ++i) {
        vector2_t intersection_parameters = intersect_ray_sphere(ray, spheres[i]);
        if (t_min <= intersection_parameters.x && intersection_parameters.x <= t_max) {
            return true;
        }

        if (t_min <= intersection_parameters.y && intersection_parameters.y <= t_max) {
            return true;
        }
    }
    return false;
}

double compute_lighting(vector3_t contact_point, vector3_t normal, vector3_t contact_point_to_camera, double specular_exponent) {
    double intensity = 0.0;
    for (size_t k = 0; k < LIGHTS_SIZE; ++k) {
        const light_t *light = lights + k;
        if (light->type == LIGHT_TYPE_AMBIENT) {
            intensity += light->intensity;
        } else {
            vector3_t light_vector = light->data.direction;
            double t_max = INFINITY;
            if (light->type == LIGHT_TYPE_POINT) {
                light_vector = vector3_subtract(light->data.position, contact_point);
                t_max = 1;
            }

            // ray_sphere_intersection_t shadow_intersection = closest_intersection(
            //     (ray_t){ .position = contact_point, .direction = light_vector },
            //     0.001,
            //     t_max);
            // if (shadow_intersection.sphere) {
            //     continue;
            // }

            if (is_shadowed((ray_t){ .position = contact_point, .direction = light_vector },
                            0.001,
                            t_max)) {
                continue;
            }

            // diffuse
            double dot = vector3_dot(light_vector, normal);
            if (dot > 0) {
                intensity += light->intensity * dot / (vector3_length(light_vector) * vector3_length(normal));
            }

            // specular
            if (specular_exponent == -1) {
                continue;
            }
            vector3_t reflection = vector3_reflect(light_vector, normal);
            double dot2 = vector3_dot(reflection, contact_point_to_camera);
            if (dot2 > 0) {
                intensity += light->intensity * pow(dot2 / (vector3_length(reflection) * vector3_length(contact_point_to_camera)), specular_exponent);
            }
        }
    }

    return intensity;
}

color_t trace_ray(ray_t ray, double t_min, double t_max, int depth) {
    ray_sphere_intersection_t intersection = closest_intersection(ray, t_min, t_max);
    if (!intersection.sphere) {
        return BACKGROUND_COLOR;
    }

    vector3_t contact_point = vector3_add(ray.position, vector3_scale(ray.direction, intersection.t));
    vector3_t normal = vector3_normalize(vector3_subtract(contact_point, intersection.sphere->center));

    double lighting = compute_lighting(contact_point, normal, vector3_negate(ray.direction), intersection.sphere->specular);
    color_t local_color = {
        .r = clamp((intersection.sphere->color.r * lighting), 0, 1),
        .g = clamp((intersection.sphere->color.g * lighting), 0, 1),
        .b = clamp((intersection.sphere->color.b * lighting), 0, 1)
    };

    double r = intersection.sphere->reflective;
    if (depth <= 0 || r <= 0) {
        return local_color;
    }

    vector3_t reflected_direction = vector3_reflect(vector3_negate(ray.direction), normal);
    color_t reflected_color = trace_ray((ray_t){ .position = contact_point, .direction = reflected_direction }, 0.001, INFINITY, depth - 1);

    return (color_t){
        .r = lerp(local_color.r, reflected_color.r, r),
        .g = lerp(local_color.g, reflected_color.g, r),
        .b = lerp(local_color.b, reflected_color.b, r)
    };
}

vector3_t canvas_to_viewport(int canvas_x, int canvas_y) {
    return (vector3_t){
        .x = (double)canvas_x * VIEWPORT_WIDTH / CANVAS_WIDTH,
        .y = (double)canvas_y * VIEWPORT_HEIGHT / CANVAS_HEIGHT,
        .z = VIEWPORT_CAMERA_DISTANCE
    };
}

vector2_t canvas_to_screen(int canvas_x, int canvas_y) {
    return (vector2_t){
        .x = CANVAS_WIDTH / 2.0 + (double)canvas_x,
        .y = CANVAS_HEIGHT / 2.0 - (double)canvas_y
    };
}

void put_pixel(int canvas_x, int canvas_y, color_t color) {
    vector2_t screen_coord = canvas_to_screen(canvas_x, canvas_y);
    DrawPixel(
        (int)screen_coord.x,
        (int)screen_coord.y,
        (Color){
            .r = (uint8_t)round(color.r * 255),
            .g = (uint8_t)round(color.g * 255),
            .b = (uint8_t)round(color.b * 255),
            .a = 255 });
}

RenderTexture target;

void draw() {
    BeginDrawing();
    Rectangle sourceRec = { 0.0, 0.0, (float)target.texture.width, -(float)target.texture.height };
    Rectangle destRec = { 0.0, 0.0, CANVAS_WIDTH, CANVAS_HEIGHT };
    DrawTexturePro(target.texture, sourceRec, destRec, (Vector2){ 0, 0 }, 0, WHITE);
    EndDrawing();
}

void update() {
    // empty
}

int main() {
    InitWindow(CANVAS_WIDTH, CANVAS_HEIGHT, "Computer Graphics from Scratch - Raytracing");

    SetTargetFPS(60);

    camera_t camera = {
        .position = { 3, 4, -2 },
        .rotation = matrix3_multiply(matrix3_y_rotation_from_angle(-M_PI / 6), matrix3_x_rotation_from_angle(M_PI / 6))
    };
    target = LoadRenderTexture(CANVAS_WIDTH, CANVAS_HEIGHT);

    BeginTextureMode(target);
    ClearBackground(BLACK);
    for (int x = -CANVAS_WIDTH / 2; x <= CANVAS_WIDTH / 2; ++x) {
        for (int y = -CANVAS_HEIGHT / 2; y <= CANVAS_HEIGHT / 2; ++y) {
            ray_t ray = { .position = camera.position,
                          .direction = vector3_transform(canvas_to_viewport(x, y), camera.rotation) };
            color_t color = trace_ray(ray, 1, INFINITY, 3);
            put_pixel(x, y, color);
        }
    }
    EndTextureMode();

    while (!WindowShouldClose()) {
        update();
        draw();
    }

    CloseWindow();

    return 0;
}
