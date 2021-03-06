/*
 * Copyright © 2012 Bart Massey
 * [This program is licensed under the "MIT License"]
 * Please see the file COPYING in the source
 * distribution of this work for license terms.
 */

/*
 * Generate Gravner-Griffeath 2D "Snowfakes"
 * 
 * Janko Gravner and David Griffeath, "Modeling Snow Crystal
 * Growth II: A mesoscopic lattice map with plausible
 * dynamics", Physica D: Nonlinear Phenomena (237)385–404
 * 2008, URL http://psoup.math.wisc.edu/papers/h2l.pdf
 * accessed 2012/12/17.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int size;
static int center;

/* initial vapor density: typ 0.3..0.9 */
static float gg_rho = 0.42;
/* freezing fraction for boundary: typ 0.001..0.02 */
static float gg_kappa = 0.01;
/* min boundary mass to join crystal for 1..2 neighbors: typ 1.05..3.0 */
static float gg_beta = 1.9;
/* max neighborhood diffusive mass to join for 3 neighbors: typ 0.01..0.04 */
static float gg_theta = 0.025;
/* min boundary mass to join for 3 neighbors: typ 0.02..0.1 */
static float gg_alpha = 0.08;
/* melting boundary mass diffusive fraction: typ "small" 0.04..0.09 */
static float gg_mu = 0.06;
/* crystal boundary mass diffusive fraction: typ "very small" */
static float gg_gamma = 0.006;
/* noise in diffusive mass: typ "tiny" */
static float gg_sigma = 0;

struct site_state {
    short attached;
    short attached_neighbors;
    float boundary_mass;
    float crystal_mass;
    float diffusive_mass;
};

static struct site_state **sites0, **sites;

static struct site_state ***sites_pages;

static int neighbor_offsets[6][2] = {
    {-1, -1},
    {-1, 0},
    {0, -1},
    {0, 1},
    {1, 0},
    {1, 1}
};

static inline int next_neighbor(int i, int r0, int c0, int *rp, int *cp) {
    int r = r0 + neighbor_offsets[i][0];
    int c = c0 + neighbor_offsets[i][1];
    if (r < 0 || r >= size || c < 0 || c >= size)
        return 0;
    *rp = r;
    *cp = c;
    return 1;
}

static void flip_sites(void) {
    struct site_state **tmp = sites;
    sites = sites0;
    sites0 = tmp;
}

static void init_sites(void) {
    struct site_state off_site = {
        .attached = 0,
        .attached_neighbors = 0,
        .boundary_mass = 0,
        .crystal_mass = 0,
        .diffusive_mass = gg_rho
    };
    struct site_state start_site = {
        .attached = 1,
        .attached_neighbors = 0,
        .boundary_mass = 0,
        .crystal_mass = 1,
        .diffusive_mass = 0
    };
    assert(size % 2 == 1);
    sites_pages = malloc(2 * sizeof(sites[0]));
    assert(sites_pages);
    for (int p = 0; p < 2; p++) {
        sites_pages[p] = malloc(size * sizeof(sites[0]));
        assert(sites_pages[p]);
        for (int r = 0; r < size; r++) {
            sites_pages[p][r] = malloc(size * sizeof(sites[r][0]));
            assert(sites_pages[p][r]);
        }
    }
    sites0 = sites_pages[0];
    sites = sites_pages[1];
    for (int r = 0; r < size; r++)
        for (int c = 0; c < size; c++)
            sites[r][c] = off_site;
    center = size / 2 + 1;
    sites[center][center] = start_site;
    for (int i = 0; i < 6; i++) {
        int r, c;
        assert(next_neighbor(i, center, center, &r, &c));
        sites[r][c].attached_neighbors = 1;
    }
    flip_sites();
}

static void diffusion(void) {
    for (int r = 1; r < size - 1; r++) {
        for (int c = 1; c < size - 1; c++) {
            /* XXX These assignments are here for convenience
               only. They should probably go somewhere else. */
            sites[r][c].attached_neighbors =
              sites0[r][c].attached_neighbors;
            sites[r][c].attached =
              sites0[r][c].attached;

            if (sites[r][c].attached) {
                sites[r][c].diffusive_mass = 0;
                continue;
            }
            sites[r][c].diffusive_mass =
              sites0[r][c].diffusive_mass / 7.0;
            for (int i = 0; i < 6; i++) {
                int rr, cc;
                if (!next_neighbor(i, r, c, &rr, &cc))
                    continue;
                /* "reflective boundary conditions": replace
                   mass from crystalline neighbors with self mass. */
                if (sites0[rr][cc].attached) {
                    rr = r;
                    cc = c;
                }
                sites[r][c].diffusive_mass +=
                  sites0[rr][cc].diffusive_mass / 7.0;
            }
        }
    }
}

static void freezing(void) {
    for (int r = 1; r < size - 1; r++) {
        for (int c = 1; c < size - 1; c++) {
            if (sites[r][c].attached || sites[r][c].attached_neighbors == 0)
                continue;
            float b0 = sites0[r][c].boundary_mass;
            float c0 = sites0[r][c].crystal_mass;
            float d0 = sites[r][c].diffusive_mass;
            sites[r][c].diffusive_mass = 0;
            sites[r][c].crystal_mass = c0 + gg_kappa * d0;
            sites[r][c].boundary_mass = b0 + (1 - gg_kappa) * d0;
        }
    }
}

static int attachment(void) {
    int stop = 0;
    for (int r = 1; r < size - 1; r++) {
        for (int c = 1; c < size - 1; c++) {
            if (sites[r][c].attached)
                continue;
            switch(sites[r][c].attached_neighbors) {
            case 0:
                /* outside boundary, so ignore */
                break;
            case 1:
            case 2:
                /* tip or edge */
#ifdef DEBUG
                printf("boundary 1/2 %d %d %f\n",
                       r, c, sites[r][c].boundary_mass);
#endif
                if (sites[r][c].boundary_mass >= gg_beta)
                    sites[r][c].attached = 1;
                break;
            case 3:
                /* concavity */
#ifdef DEBUG
                printf("boundary 3 %d %d %f %f\n",
                       r, c,
                       sites[r][c].boundary_mass,
                       sites[r][c].diffusive_mass);
#endif
                if (sites[r][c].boundary_mass >= 1.0) {
                    sites[r][c].attached = 1;
                    break;
                }
                if (sites[r][c].boundary_mass < gg_alpha)
                    break;
                float diffusive_mass = sites[r][c].diffusive_mass;
                for (int i = 0; i < 6; i++) {
                    int rr, cc;
                    if (!next_neighbor(i, r, c, &rr, &cc))
                        continue;
                    diffusive_mass += sites[rr][cc].diffusive_mass;
                }
                if (diffusive_mass < gg_theta)
                    sites[r][c].attached = 1;
                break;
            default:
                /* hole, essentially */
                sites[r][c].attached = 1;
            }
            if (sites[r][c].attached) {
#ifdef DEBUG
                printf("attached %d %d\n", r, c);
#endif
                for (int i = 0; i < 6; i++) {
                    int rr, cc;
                    if (!next_neighbor(i, r, c, &rr, &cc))
                        continue;
                    sites[rr][cc].attached_neighbors++;
                }
                sites[r][c].crystal_mass += sites[r][c].boundary_mass;
                sites[r][c].boundary_mass = 0;
                int ss = size / 3;
                if (r < ss || r >= size - ss || c < ss || c >= size - ss)
                    stop = 1;
            }
        }
    }
    return stop;
}

static void melting(void) {
    for (int r = 1; r < size - 1; r++) {
        for (int c = 1; c < size - 1; c++) {
            if (sites[r][c].attached || sites[r][c].attached_neighbors == 0)
                continue;
            float b0 = sites[r][c].boundary_mass;
            float c0 = sites[r][c].crystal_mass;
            float d0 = sites[r][c].diffusive_mass;
            sites[r][c].boundary_mass = (1 - gg_mu) * b0;
            sites[r][c].crystal_mass = (1 - gg_gamma) * c0;
            sites[r][c].diffusive_mass = d0 + gg_mu * b0 + gg_gamma * c0;
        }
    }
}

static void noise(void) {
    for (int r = 1; r < size - 1; r++) {
        for (int c = 1; c < size - 1; c++) {
            double sigma = gg_sigma;
            if ((random() & 1) == 0)
                sigma = -sigma;
            sites[r][c].diffusive_mass *= (1 + sigma);
        }
    }
}

static void render(void) {
    float scale = 1000.0;
    float yscale = 1.0 / sqrtf(3.0);
    float dscale = scale / size;
    float dotscale = 0.25 * dscale;
    printf("<?xml version=\"1.0\"?>\n"
           "<svg width=\"%f\" height=\"%f\"\n"
           "version=\"1.1\"\n"
           "xmlns=\"http://www.w3.org/2000/svg\">\n",
           scale, scale * yscale);
    for (int r = 0; r < size; r++) {
        for (int c = 0; c < size; c++) {
            if (!sites[r][c].attached)
                continue;
            int x0 = r - center;
            int y0 = c - center;
            float d = hypotf(x0, y0);
            float a = atan2f(y0, x0);
            a += M_PI_4;
            float x = (d * cosf(a) + center) * dscale;
            float y = (d * sinf(a) + center) * dscale * yscale;
            printf("  <circle cx=\"%f\" cy=\"%f\" r=\"%f\"/>\n",
                   x, y, sites[r][c].crystal_mass * dotscale);
        }
    }
    printf("</svg>\n");
}

int main(int argc, char **argv) {
    assert(argc == 2);
    size = atoi(argv[1]);
    assert(size > 0);
    init_sites();
    int t = 1;
    while (1) {
        if (t % 1000 == 0)
            fprintf(stderr, ".");
        if (t >= 100000) {
            fprintf(stderr, "!");
            break;
        }
        diffusion();
        freezing();
        int stop = attachment();
        if (stop)
            break;
        melting();
        if (gg_gamma > 0)
            noise();
        flip_sites();
        t++;
    }
    fprintf(stderr, "\n");
    render();
    return 0;
}
