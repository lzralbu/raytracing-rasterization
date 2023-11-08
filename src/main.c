#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

typedef struct Sphere {
    Vector3 center;
    float radius;
    Color color;
    float specular;
} Sphere;

enum LIGHT_TYPE {
    LIGHT_TYPE_POINT,
    LIGHT_TYPE_DIRECTIONAL,
    LIGHT_TYPE_AMBIENT
};

typedef struct Light {
    int type;
    float intensity;
    union {
        Vector3 position;
        Vector3 direction;
    } data;
} Light;

static const int CANVAS_WIDTH = 800;
static const int CANVAS_HEIGHT = 800;

static const float VIEWPORT_WIDTH = 1;
static const float VIEWPORT_HEIGHT = 1;
static const float VIEWPORT_CAMERA_DISTANCE = 1;

#define BACKGROUND_COLOR RAYWHITE

static Sphere spheres[] = {
    [0] = { .center = { 0, -1, 3 }, .radius = 1, .color = { 255, 0, 0, 255 }, .specular = 500 },
    [1] = { .center = { 2, 0, 4 }, .radius = 1, .color = { 0, 0, 255, 255 }, .specular = 500 },
    [2] = { .center = { -2, 0, 4 }, .radius = 1, .color = { 0, 255, 0, 255 }, .specular = 10 },
    [3] = { .center = { 0, -5001, 0 }, .radius = 5000, .color = { 255, 255, 0, 255 }, .specular = 1000 }
};
#define SPHERES_SIZE (sizeof(spheres) / sizeof(Sphere))

static Light lights[] = {
    [0] = { .type = LIGHT_TYPE_AMBIENT, .intensity = 0.2f },
    [1] = { .type = LIGHT_TYPE_POINT, .intensity = 0.6f, .data.position = { 2, 1, 0 } },
    [2] = { .type = LIGHT_TYPE_DIRECTIONAL, .intensity = 0.2f, .data.direction = { 1, 4, 4 } }
};
#define LIGHTS_SIZE (sizeof(lights) / sizeof(Light))

Vector2 intersect_ray_sphere(Ray ray, Sphere sphere) {
    Vector3 temp_vec = Vector3Subtract(ray.position, sphere.center);

    float a = Vector3DotProduct(ray.direction, ray.direction);
    float b = 2 * Vector3DotProduct(temp_vec, ray.direction);
    float c = Vector3DotProduct(temp_vec, temp_vec) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return (Vector2){ .x = INFINITY, .y = INFINITY };
    }

    float t1 = (-b + sqrtf(discriminant)) / (2 * a);
    float t2 = (-b - sqrtf(discriminant)) / (2 * a);
    return (Vector2){ .x = t1, .y = t2 };
}

float compute_lighting(Vector3 contact_point, Vector3 normal, Vector3 contact_point_to_camera, float specular_exponent) {
    float intensity = 0.0f;
    for (size_t k = 0; k < LIGHTS_SIZE; ++k) {
        const Light *light = lights + k;
        if (light->type == LIGHT_TYPE_AMBIENT) {
            intensity += light->intensity;
        } else {
            Vector3 light_vector =
                light->type == LIGHT_TYPE_POINT
                    ? Vector3Subtract(light->data.position, contact_point)
                    : light->data.direction;

            // diffuse
            float dot = Vector3DotProduct(light_vector, normal);
            if (dot > 0) {
                intensity += light->intensity * dot / (Vector3Length(light_vector) * Vector3Length(normal));
            }

            // specular
            if (specular_exponent == -1) {
                continue;
            }
            Vector3 reflection = Vector3Subtract(Vector3Scale(normal, 2 * Vector3DotProduct(normal, light_vector)), light_vector);
            float dot2 = Vector3DotProduct(reflection, contact_point_to_camera);
            if (dot2 > 0) {
                intensity += light->intensity * powf(dot2 / (Vector3Length(reflection) * Vector3Length(contact_point_to_camera)), specular_exponent);
            }
        }
    }

    return intensity;
}

Color trace_ray(Ray ray, float t_min, float t_max) {
    float closest_t = INFINITY;
    Sphere const *closest_sphere = 0;
    for (size_t i = 0; i < SPHERES_SIZE; ++i) {
        Vector2 intersection_parameters = intersect_ray_sphere(ray, spheres[i]);
        if (t_min <= intersection_parameters.x && intersection_parameters.x <= t_max && intersection_parameters.x < closest_t) {
            closest_t = intersection_parameters.x;
            closest_sphere = spheres + i;
        }

        if (t_min <= intersection_parameters.y && intersection_parameters.y <= t_max && intersection_parameters.y < closest_t) {
            closest_t = intersection_parameters.y;
            closest_sphere = spheres + i;
        }
    }
    if (closest_sphere == 0) {
        return BACKGROUND_COLOR;
    }

    Vector3 contact_point = Vector3Add(ray.position, Vector3Scale(ray.direction, closest_t));
    Vector3 normal = Vector3Normalize(Vector3Subtract(contact_point, closest_sphere->center));

    float lighting = compute_lighting(contact_point, normal, Vector3Negate(ray.direction), closest_sphere->specular);
    Color new_color = {
        .r = (uint8_t)Clamp((closest_sphere->color.r * lighting), 0, 255),
        .g = (uint8_t)Clamp((closest_sphere->color.g * lighting), 0, 255),
        .b = (uint8_t)Clamp((closest_sphere->color.b * lighting), 0, 255),
        .a = 255
    };

    return new_color;
}

void put_pixel(int canvas_x, int canvas_y, Color color) {
    int screen_x = CANVAS_WIDTH / 2 + canvas_x;
    int screen_y = CANVAS_HEIGHT / 2 - canvas_y;
    DrawPixel(screen_x, screen_y, color);
}

Vector3 canvas_to_viewport(int canvas_x, int canvas_y) {
    return (Vector3){
        .x = (float)canvas_x * VIEWPORT_WIDTH / CANVAS_WIDTH,
        .y = (float)canvas_y * VIEWPORT_HEIGHT / CANVAS_HEIGHT,
        .z = VIEWPORT_CAMERA_DISTANCE
    };
}

void draw() {
    BeginDrawing();
    ClearBackground(BACKGROUND_COLOR);

    // rlPushMatrix();
    // rlMatrixMode(RL_PROJECTION);
    // rlLoadIdentity();

    // rlTranslatef(CANVAS_WIDTH / 2, CANVAS_HEIGHT / 2, 0);
    // rlScalef(1, -1, 1);

    Ray ray = { 0 };
    for (int x = -CANVAS_WIDTH / 2; x <= CANVAS_WIDTH / 2; ++x) {
        for (int y = -CANVAS_HEIGHT / 2; y <= CANVAS_HEIGHT / 2; ++y) {
            ray.direction = canvas_to_viewport(x, y);
            Color color = trace_ray(ray, 1, INFINITY);
            put_pixel(x, y, color);
            // DrawPixel(x, y, color);
        }
    }

    // rlPopMatrix();

    EndDrawing();
}

void update() {
    // empty
}

int main() {
    InitWindow(CANVAS_WIDTH, CANVAS_HEIGHT, "Computer Graphics from Scratch - Raytracing");

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        update();
        draw();
    }

    CloseWindow();

    return 0;
}
