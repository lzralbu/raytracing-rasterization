#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

double clamp(double value, double lower, double upper) {
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

typedef struct vector2_t {
    double x;
    double y;
} vector2_t;

typedef struct vector3_t {
    double x;
    double y;
    double z;
} vector3_t;

typedef struct color_t {
    double r;
    double g;
    double b;
} color_t;

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
    [0] = { .type = LIGHT_TYPE_AMBIENT, .intensity = 0.2f },
    [1] = { .type = LIGHT_TYPE_POINT, .intensity = 0.6f, .data.position = { 2, 1, 0 } },
    [2] = { .type = LIGHT_TYPE_DIRECTIONAL, .intensity = 0.2f, .data.direction = { 1, 4, 4 } }
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

            ray_sphere_intersection_t shadow_intersection = closest_intersection(
                (ray_t){ .position = contact_point, .direction = light_vector },
                0.001,
                t_max);
            if (shadow_intersection.sphere) {
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
        .r = (1.0 - r) * local_color.r + r * reflected_color.r,
        .g = (1.0 - r) * local_color.g + r * reflected_color.g,
        .b = (1.0 - r) * local_color.b + r * reflected_color.b,
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

    target = LoadRenderTexture(CANVAS_WIDTH, CANVAS_HEIGHT);

    BeginTextureMode(target);
    ClearBackground(BLACK);
    ray_t ray = { 0 };
    for (int x = -CANVAS_WIDTH / 2; x <= CANVAS_WIDTH / 2; ++x) {
        for (int y = -CANVAS_HEIGHT / 2; y <= CANVAS_HEIGHT / 2; ++y) {
            ray.direction = canvas_to_viewport(x, y);
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
