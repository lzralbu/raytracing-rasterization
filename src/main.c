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
} Sphere;

const int CANVAS_WIDTH = 600;
const int CANVAS_HEIGHT = 600;

const float VIEWPORT_WIDTH = 1;
const float VIEWPORT_HEIGHT = 1;
const float VIEWPORT_CAMERA_DISTANCE = 1;

#define BACKGROUND_COLOR RAYWHITE

#define SPHERES_SIZE (size_t)3
Sphere spheres[SPHERES_SIZE] = {
    [0] = { .center = { 0, -1, 3 }, .radius = 1, .color = { 255, 0, 0, 255 } },
    [1] = { .center = { 2, 0, 4 }, .radius = 1, .color = { 0, 0, 255, 255 } },
    [2] = { .center = { -2, 0, 4 }, .radius = 1, .color = { 0, 255, 0, 255 } },
};

Vector3 canvas_to_viewport(int canvas_x, int canvas_y) {
    return (Vector3){
        .x = (float)canvas_x * VIEWPORT_WIDTH / CANVAS_WIDTH,
        .y = (float)canvas_y * VIEWPORT_HEIGHT / CANVAS_HEIGHT,
        .z = VIEWPORT_CAMERA_DISTANCE
    };
}

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

    return closest_sphere->color;
}

void put_pixel(int canvas_x, int canvas_y, Color color) {
    int screen_x = CANVAS_WIDTH / 2 + canvas_x;
    int screen_y = CANVAS_HEIGHT / 2 - canvas_y;
    DrawPixel(screen_x, screen_y, color);
}

void draw() {
    BeginDrawing();
    ClearBackground(BACKGROUND_COLOR);

    Ray ray = { 0 };
    for (int x = -CANVAS_WIDTH / 2; x <= CANVAS_WIDTH / 2; ++x) {
        for (int y = -CANVAS_HEIGHT / 2; y <= CANVAS_HEIGHT / 2; ++y) {
            ray.direction = canvas_to_viewport(x, y);
            Color color = trace_ray(ray, 1, INFINITY);
            put_pixel(x, y, color);
        }
    }

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
