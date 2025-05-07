#include "devoir_3.h"

#include <stdio.h>
#include <stdlib.h>


int read_initial_conditions(const char *filename, double **u, double **v) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", filename);
        return -1;
    }
    
    char line[256];
    int n = 0;
    while (fgets(line, sizeof(line), file) != NULL) n++;
    rewind(file);
    
    double *_u = (double *) malloc(2 * n * sizeof(double));
    double *_v = (double *) malloc(2 * n * sizeof(double)); 
    if (_u == NULL || _v == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }

    int read;
    int i = 0;
    while ((read = fscanf(file, "%lf %lf %lf %lf", &_u[2*i], &_u[2*i+1], &_v[2*i], &_v[2*i+1])) != EOF) {
        if (read != 4) {
            fprintf(stderr, "Error: Incorrect format in file %s\n", filename);
            fprintf(stderr, "Expected '<ux_i> <uy_i> <vx_i> <vy_i>'\n");
            fclose(file);
            free(_u);
            free(_v);
            return -1;
        }

        i++;
    }

    if (i != n) {
        fprintf(stderr, "Error: Number of lines read does not match expected count\n");
        fclose(file);
        free(_u);
        free(_v);
        return -1;
    }

    *u = _u;
    *v = _v;

    fclose(file);
    
    return n;
}

int newmark(
    const SymBandMatrix *Kbd, const SymBandMatrix *Mbd, double *u, double *v,
    const int n, const double dt, const double T, const double gamma, const double beta        
) {

    CSRMatrix *M_csr = band_to_csr(Mbd);
    CSRMatrix *K_csr = band_to_csr(Kbd);

    // assemble_system guarantees that both matrices are the same size and same band.
    SymBandMatrix *lhs_bd = allocate_sym_band_matrix(Kbd->n, Kbd->k);
    for (int k = 0; k < lhs_bd->n * (lhs_bd->k + 1); k++)
        lhs_bd->data[k] = Mbd->data[k] + dt*dt*beta * Kbd->data[k];

    // Factorize the left-hand side matrix
    sym_band_LDL(lhs_bd->data, lhs_bd->n, lhs_bd->k);

    SymBandMatrix *rhs_bd = allocate_sym_band_matrix(Kbd->n, Kbd->k);
    for (int k = 0; k < rhs_bd->n * (rhs_bd->k + 1); k++)
        rhs_bd->data[k] = Mbd->data[k] - (dt*dt/2.)*(1.-2.*beta) * Kbd->data[k];

    // For future matrix multiplications
    CSRMatrix *rhs_csr = band_to_csr(rhs_bd);

    double *p = (double *) malloc(2 * n * sizeof(double));
    double *rhs = (double *) malloc(2 * n * sizeof(double));
    double *interp = (double *) malloc(2 * n * sizeof(double));
    if (p == NULL || rhs == NULL || interp == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free_csr(M_csr);
        free_csr(K_csr);
        free_csr(rhs_csr);
        free_sym_band_matrix(lhs_bd);
        free_sym_band_matrix(rhs_bd);
        return -1;
    }
    
    // p = Mv
    Matvec(M_csr->n, M_csr->row_ptr, M_csr->col_idx, M_csr->data, v, p);

    double t = 0.;
    while (t < T) {
        printf("%lf\n", u[0]);
        if (fabs(u[0]) > 1e12) {
            fprintf(stderr, "Fatal error: Displacement exceeds threshold\n");
            exit(EXIT_FAILURE);
        }

        // Compute the right-hand side of first equation
        Matvec(rhs_csr->n, rhs_csr->row_ptr, rhs_csr->col_idx, rhs_csr->data, u, rhs);
        for (int k = 0; k < 2 * n; k++)
            rhs[k] += dt * p[k];

        // Solve the system, new u is in rhs 
        solve_sym_band(lhs_bd->data, lhs_bd->n, lhs_bd->k, rhs);

        // Linear interpolation of u
        for (int k = 0; k < 2 * n; k++)
            interp[k] = (1-gamma) * u[k] + gamma * rhs[k];

        // here, v just serve as a temporary variable
        Matvec(K_csr->n, K_csr->row_ptr, K_csr->col_idx, K_csr->data, interp, v);
        
        // Update velocity
        for (int k = 0; k < 2 * n; k++)
            p[k] -= dt * v[k]; 
        
        memcpy(u, rhs, 2 * n * sizeof(double));

        t += dt;
    }

    SymBandMatrix Mbd_cpy;
    memcpy(&Mbd_cpy, Mbd, sizeof(SymBandMatrix));
    sym_band_LDL(Mbd_cpy.data, Mbd_cpy.n, Mbd_cpy.k);
    solve_sym_band(Mbd_cpy.data, Mbd_cpy.n, Mbd_cpy.k, p);
    memcpy(v, p, 2 * n * sizeof(double));

    free_csr(M_csr);
    free_csr(K_csr);
    free_csr(rhs_csr);
    free_sym_band_matrix(lhs_bd);
    free_sym_band_matrix(rhs_bd);
    free(p);
    free(rhs);
    free(interp);

    return 0;
}
