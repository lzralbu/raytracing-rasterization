#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "raylib.h"
#include "vec.h"

#include "vmath.h"

enum GEOMETRY_TYPE { GEOMETRY_TYPE_SPHERE, GEOMETRY_TYPE_TRIANGLE };

typedef struct geometry_t {
    int type;
    void *data;
} geometry_t;

typedef struct geometry_data_sphere_t {
    vector3_t center;
    double radius;
} geometry_data_sphere_t;

typedef struct geometry_data_triangle_t {
    vector3_t vertices[3];
} geometry_data_triangle_t;

void geometry_get_normal_at(vector3_t *result, geometry_t const *geometry, vector3_t const *point) {
    if (geometry->type == GEOMETRY_TYPE_SPHERE) {
        geometry_data_sphere_t const *sphere = (geometry_data_sphere_t const *)geometry->data;
        vector3_t diff;
        vector3_subtract(&diff, point, &sphere->center);
        vector3_normalize(result, &diff);
    } else if (geometry->type == GEOMETRY_TYPE_TRIANGLE) {
        geometry_data_triangle_t const *triangle = (geometry_data_triangle_t const *)geometry->data;
        vector3_t ab;
        vector3_subtract(&ab, &triangle->vertices[1], &triangle->vertices[0]);
        vector3_t ac;
        vector3_subtract(&ac, &triangle->vertices[2], &triangle->vertices[0]);
        vector3_t cross;
        vector3_cross(&cross, &ab, &ac);
        vector3_normalize(result, &cross);
    } else {
        result->a[0] = INFINITY;
        result->a[1] = INFINITY;
        result->a[2] = INFINITY;
    }
}

typedef struct material_t {
    color_t color;
    double specular;
    double reflective;
} material_t;

typedef struct primitive_t {
    geometry_t geometry;
    material_t material;
} primitive_t;

typedef struct ray_t {
    vector3_t position;
    vector3_t direction;
} ray_t;

typedef struct ray_primitive_intersection_t {
    primitive_t *primitive;
    vec_double_t parameters;
} ray_primitive_intersection_t;

enum LIGHT_TYPE { LIGHT_TYPE_POINT, LIGHT_TYPE_DIRECTIONAL, LIGHT_TYPE_AMBIENT };

typedef struct light_t {
    int type;
    double intensity;

    union {
        vector3_t position;
        vector3_t direction;
    } data;
} light_t;

typedef struct camera_t {
    vector3_t position;
    matrix_t rotation;
} camera_t;

typedef vec_t(primitive_t) primitive_vec_t;
typedef vec_t(light_t) light_vec_t;

typedef struct scene_t {
    camera_t camera;
    primitive_vec_t primitives;
    light_vec_t lights;
} scene_t;

void scene_init(scene_t *scene) {
    scene->camera.position = (vector3_t){ 0, 2, -2 };

    matrix_t y_rotation;
    matrix_y_rotation_from_angle(&y_rotation, 0);
    matrix_t x_rotation;
    matrix_x_rotation_from_angle(&x_rotation, M_PI / 6);
    matrix_multiply(&scene->camera.rotation, &y_rotation, &x_rotation);

    static geometry_data_sphere_t geometry_data_spheres[] = { [0] = {
                                                         .center = { 0, -1, 3 },
                                                         .radius = 1,
                                                     },
                                                     [1] = {
                                                        .center = { 2, 0, 4 },
                                                        .radius = 1,
                                                     },
                                                     [2] = {
                                                        .center = { -2, 0, 4 },
                                                        .radius = 1,
                                                     },
                                                     [3] = {
                                                        .center = { 0, -5001, 0 },
                                                        .radius = 5000,
                                                     } };

    static geometry_data_triangle_t geometry_data_triangles[] = {
        [0] = { .vertices[0] = { .a[0] = 0, .a[1] = 1, .a[2] = 6 },
                .vertices[1] = { .a[0] = -0.5, .a[1] = 0, .a[2] = 6 },
                .vertices[2] = { .a[0] = 0.5, .a[1] = 0, .a[2] = 6 } },
        [1] = { .vertices[0] = { .a[0] = -1, .a[1] = 1, .a[2] = 1 },
                .vertices[1] = { .a[0] = -1, .a[1] = 0, .a[2] = 3 },
                .vertices[2] = { .a[0] = 0, .a[1] = 0, .a[2] = 1 } }
    };

    material_t materials[] = { [0] = { .color = { 1, 0, 0, 1 }, .specular = 500, .reflective = 0.2 },
                               [1] = { .color = { 0, 0, 1, 1 }, .specular = 500, .reflective = 0.3 },
                               [2] = { .color = { 0, 1, 0, 1 }, .specular = 10, .reflective = 0.4 },
                               [3] = { .color = { 1, 1, 0, 1 }, .specular = 1000, .reflective = 0.5 } };

    primitive_t primitives[] = {
        [0] = { .geometry = { .type = GEOMETRY_TYPE_SPHERE, .data = &geometry_data_spheres[0] },
                .material = materials[0] },
        [1] = { .geometry = { .type = GEOMETRY_TYPE_SPHERE, .data = &geometry_data_spheres[1] },
                .material = materials[1] },
        [2] = { .geometry = { .type = GEOMETRY_TYPE_SPHERE, .data = &geometry_data_spheres[2] },
                .material = materials[2] },
        [3] = { .geometry = { .type = GEOMETRY_TYPE_SPHERE, .data = &geometry_data_spheres[3] },
                .material = materials[3] },
        [4] = { .geometry = { .type = GEOMETRY_TYPE_TRIANGLE, .data = &geometry_data_triangles[0] },
                .material = materials[1] },
        [5] = { .geometry = { .type = GEOMETRY_TYPE_TRIANGLE, .data = &geometry_data_triangles[1] },
                .material = materials[0] }
    };
    size_t primitives_size = sizeof(primitives) / sizeof(primitive_t);

    vec_init(&scene->primitives);
    vec_reserve(&scene->primitives, primitives_size);
    for (size_t i = 0; i < primitives_size; ++i) {
        vec_push(&scene->primitives, primitives[i]);
    }

    static light_t lights[] = { [0] = { .type = LIGHT_TYPE_AMBIENT, .intensity = 0.2 },
                                [1] = { .type = LIGHT_TYPE_POINT, .intensity = 0.6, .data.position = { 2, 1, 0 } },
                                [2] = {
                                    .type = LIGHT_TYPE_DIRECTIONAL, .intensity = 0.2, .data.direction = { 1, 4, 4 } } };
    const size_t lights_size = sizeof(lights) / sizeof(light_t);

    vec_init(&scene->lights);
    vec_reserve(&scene->lights, lights_size);
    for (size_t i = 0; i < lights_size; ++i) {
        vec_push(&scene->lights, lights[i]);
    }
}

void scene_quit(scene_t *scene) {
    vec_deinit(&scene->primitives);
    vec_deinit(&scene->lights);
}

static const int CANVAS_WIDTH = 900;
static const int CANVAS_HEIGHT = 900;

static const double VIEWPORT_WIDTH = 1;
static const double VIEWPORT_HEIGHT = 1;
static const double VIEWPORT_CAMERA_DISTANCE = 1;

const color_t BACKGROUND_COLOR = { 0, 0, 0, 1 };

// const color_t BACKGROUND_COLOR = { 1, 1, 1, 1 };

static scene_t main_scene;

typedef double matrix_3x4_t[3][4];

void print_matrix(matrix_3x4_t a) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            printf("%g ", a[i][j]);
        }
        puts("");
    }
}

void swap_rows(matrix_3x4_t a, size_t row, size_t other_row) {
    if (row == other_row) {
        return;
    }

    double temp[4];
    memcpy(temp, a[row], 4 * sizeof(double));
    memcpy(a[row], a[other_row], 4 * sizeof(double));
    memcpy(a[other_row], temp, 4 * sizeof(double));
}

void gaussian_elimination(matrix_3x4_t a, size_t rows, double *solution) {
    size_t cols = rows + 1;
    for (size_t k = 1; k < rows; ++k) {
        if (a[k - 1][k - 1] == 0) {
            for (size_t i = k; i < rows; ++i) {
                if (a[i][i] != 0) {
                    swap_rows(a, i, k - 1);
                    break;
                }
            }
        }
        for (size_t i = k; i < rows; ++i) {
            double factor = a[i][k - 1] / a[k - 1][k - 1];
            for (size_t j = 0; j < cols; ++j) {
                a[i][j] -= factor * a[k - 1][j];
            }
        }

        // printf("-----------------------\n%zu step\n", k);
        // print_matrix(a);
    }

    for (size_t i = 0; i < rows; ++i) {
        solution[rows - i - 1] = a[rows - i - 1][cols - 1];
        for (size_t j = rows - i; j + 1 < cols; ++j) {
            solution[rows - i - 1] -= a[rows - i - 1][j] * solution[j];
        }
        solution[rows - i - 1] /= a[rows - i - 1][rows - i - 1];
    }
}

bool test_gaussian_elimination() {
    double a[3][4] = {
        [0] = { 1, 3, 1, 1 },
        [1] = { 2, 6, 9, 7 },
        [2] = { 2, 8, 8, 6 },
    };

    // print_matrix(a);
    // puts("-----------------");
    // swap_rows(a, 2, 0);
    // print_matrix(a);

    double solution[3];

    gaussian_elimination(a, 3, solution);

    printf("solution = (%g, %g, %g)\n", solution[0], solution[1], solution[2]);

    return true;
}

vec_double_t intersect_ray_primitive(ray_t const *ray, primitive_t const *primitive) {
    vec_double_t parameters;
    vec_init(&parameters);

    if (primitive->geometry.type == GEOMETRY_TYPE_SPHERE) {
        geometry_data_sphere_t const *geometry_data = (geometry_data_sphere_t const *)primitive->geometry.data;
        vector3_t temp_vec;
        vector3_subtract(&temp_vec, &ray->position, &geometry_data->center);
        double a = vector3_dot(&ray->direction, &ray->direction);
        double b = 2 * vector3_dot(&temp_vec, &ray->direction);
        double c = vector3_dot(&temp_vec, &temp_vec) - geometry_data->radius * geometry_data->radius;
        double discriminant = b * b - 4 * a * c;
        if (discriminant > 0) {
            vec_reserve(&parameters, 2);
            vec_push(&parameters, (-b + sqrt(discriminant)) / (2.0 * a));
            vec_push(&parameters, (-b - sqrt(discriminant)) / (2.0 * a));
        }
    } else if (primitive->geometry.type == GEOMETRY_TYPE_TRIANGLE) {
        geometry_data_triangle_t const *geometry_data = (geometry_data_triangle_t const *)primitive->geometry.data;
        vector3_t p;
        vector3_subtract(&p, &geometry_data->vertices[1], &geometry_data->vertices[0]);
        vector3_t q;
        vector3_subtract(&q, &geometry_data->vertices[2], &geometry_data->vertices[0]);
        vector3_t cross_pq;
        vector3_cross(&cross_pq, &p, &q);
        if (fabs(vector3_dot(&cross_pq, &ray->direction)) > 0.001) {
            vector3_t w;
            vector3_subtract(&w, &ray->position, &geometry_data->vertices[0]);
            double augmented_matrix[3][4] = {
                [0] = { [0] = p.a[0], [1] = q.a[0], [2] = -ray->direction.a[0], [3] = w.a[0] },
                [1] = { [0] = p.a[1], [1] = q.a[1], [2] = -ray->direction.a[1], [3] = w.a[1] },
                [2] = { [0] = p.a[2], [1] = q.a[2], [2] = -ray->direction.a[2], [3] = w.a[2] }
            };

            double solution[3];
            gaussian_elimination(augmented_matrix, 3, solution);

            if (solution[0] >= 0 && solution[1] >= 0 && solution[0] + solution[1] <= 1) {
                vec_push(&parameters, solution[2]);
            }
        }
    }

    return parameters;
}

ray_primitive_intersection_t closest_intersection(ray_t const *ray, double t_min, double t_max) {
    ray_primitive_intersection_t closest_intersection = { 0 };
    closest_intersection.primitive = 0;
    vec_init(&closest_intersection.parameters);
    vec_push(&closest_intersection.parameters, INFINITY);

    for (int i = 0; i < main_scene.primitives.length; ++i) {
        vec_double_t parameters = intersect_ray_primitive(ray, &main_scene.primitives.data[i]);
        for (int j = 0; j < parameters.length; ++j) {
            if (t_min <= parameters.data[j] && parameters.data[j] <= t_max &&
                parameters.data[j] < *closest_intersection.parameters.data) {
                *closest_intersection.parameters.data = parameters.data[j];
                closest_intersection.primitive = &main_scene.primitives.data[i];
            }
        }
        vec_deinit(&parameters);
    }

    return closest_intersection;
}

bool is_shadowed(ray_t const *ray, double t_min, double t_max) {
    for (int i = 0; i < main_scene.primitives.length; ++i) {
        vec_double_t parameters = intersect_ray_primitive(ray, &main_scene.primitives.data[i]);
        for (int j = 0; j < parameters.length; ++j) {
            if (t_min <= parameters.data[j] && parameters.data[j] <= t_max) {
                vec_deinit(&parameters);
                return true;
            }
        }
        vec_deinit(&parameters);
    }
    return false;
}

double compute_lighting(
    vector3_t const *contact_point, vector3_t const *normal, vector3_t const *contact_point_to_camera,
    double specular_exponent
) {
    double intensity = 0.0;
    for (int k = 0; k < main_scene.lights.length; ++k) {
        light_t const *light = &main_scene.lights.data[k];
        if (light->type == LIGHT_TYPE_AMBIENT) {
            intensity += light->intensity;
        } else {
            vector3_t light_vector = light->data.direction;
            double t_max = INFINITY;
            if (light->type == LIGHT_TYPE_POINT) {
                vector3_subtract(&light_vector, &light->data.position, contact_point);
                t_max = 1;
            }

            ray_t temp_ray = { .position = *contact_point, .direction = light_vector };
            if (is_shadowed(&temp_ray, 0.001, t_max)) {
                continue;
            }

            // diffuse

            double dot = vector3_dot(&light_vector, normal);
            if (dot > 0) {
                intensity += light->intensity * dot / vector3_length(&light_vector);
            }

            // specular
            if (specular_exponent == -1) {
                continue;
            }
            vector3_t reflection;
            vector3_reflect(&reflection, &light_vector, normal);
            double dot2 = vector3_dot(&reflection, contact_point_to_camera);
            if (dot2 > 0) {
                intensity += light->intensity *
                             pow(dot2 / (vector3_length(&reflection) * vector3_length(contact_point_to_camera)),
                                 specular_exponent);
            }
        }
    }

    return intensity;
}

void trace_ray(color_t *result, ray_t const *ray, double t_min, double t_max, int depth) {
    ray_primitive_intersection_t intersection = closest_intersection(ray, t_min, t_max);
    if (!intersection.primitive) {
        memcpy(result, &BACKGROUND_COLOR, sizeof(color_t));
        return;
    }

    // lights

    vector3_t ray_direction_scaled;
    vector3_scale(&ray_direction_scaled, &ray->direction, intersection.parameters.data[0]);
    vec_deinit(&intersection.parameters);

    vector3_t contact_point;
    vector3_add(&contact_point, &ray->position, &ray_direction_scaled);
    vector3_t normal;
    geometry_get_normal_at(&normal, &intersection.primitive->geometry, &contact_point);
    vector3_t ray_direction_negated;
    vector3_negate(&ray_direction_negated, &ray->direction);
    double lighting =
        compute_lighting(&contact_point, &normal, &ray_direction_negated, intersection.primitive->material.specular);
    color_t material_color_scaled;
    vector4_scale(&material_color_scaled, &intersection.primitive->material.color, lighting);
    color_t local_color;
    vector4_clamp(&local_color, &material_color_scaled, 0, 1);

    // reflection

    double r = intersection.primitive->material.reflective;
    if (depth <= 0 || r <= 0) {
        memcpy(result, &local_color, sizeof(color_t));
        return;
    }

    vector3_t reflected_direction;
    vector3_reflect(&reflected_direction, &ray_direction_negated, &normal);
    ray_t new_ray = { .position = contact_point, .direction = reflected_direction };
    color_t reflected_color;
    trace_ray(&reflected_color, &new_ray, 0.001, INFINITY, depth - 1);

    vector4_lerp(result, &local_color, &reflected_color, r);
}

void canvas_to_viewport(vector3_t *result, int canvas_x, int canvas_y) {
    result->a[0] = (double)canvas_x * VIEWPORT_WIDTH / CANVAS_WIDTH;
    result->a[1] = (double)canvas_y * VIEWPORT_HEIGHT / CANVAS_HEIGHT;
    result->a[2] = VIEWPORT_CAMERA_DISTANCE;
}

void canvas_to_screen(vector2_t *result, int canvas_x, int canvas_y) {
    result->a[0] = CANVAS_WIDTH / 2.0 + (double)canvas_x;
    result->a[1] = CANVAS_HEIGHT / 2.0 - (double)canvas_y;
}

void put_pixel(int canvas_x, int canvas_y, color_t color) {
    vector2_t screen_coord;
    canvas_to_screen(&screen_coord, canvas_x, canvas_y);

    DrawPixel(
        (int)screen_coord.a[0],
        (int)screen_coord.a[1],
        (Color){ .r = (uint8_t)round(color.a[0] * 255),
                 .g = (uint8_t)round(color.a[1] * 255),
                 .b = (uint8_t)round(color.a[2] * 255),
                 .a = 255 }
    );
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

    target = LoadRenderTexture(CANVAS_WIDTH, CANVAS_HEIGHT);

    test_gaussian_elimination();

    scene_init(&main_scene);

    BeginTextureMode(target);
    ClearBackground(BLACK);
    for (int x = -CANVAS_WIDTH / 2; x <= CANVAS_WIDTH / 2; ++x) {
        for (int y = -CANVAS_HEIGHT / 2; y <= CANVAS_HEIGHT / 2; ++y) {
            vector3_t viewport_coords;
            canvas_to_viewport(&viewport_coords, x, y);

            ray_t ray;
            ray.position = main_scene.camera.position;
            vector3_transform(&ray.direction, &viewport_coords, &main_scene.camera.rotation);

            color_t color;
            trace_ray(&color, &ray, 1, INFINITY, 3);

            put_pixel(x, y, color);
        }
    }
    EndTextureMode();

    scene_quit(&main_scene);

    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update();
        draw();
    }

    CloseWindow();

    return 0;
}
