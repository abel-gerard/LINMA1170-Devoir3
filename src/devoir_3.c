#include "devoir_3.h"


#define DEV_3_LOG 1

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
    const unsigned int n, const double dt, const double T, const double gamma, const double beta,
    const char *final_filename, const char *time_filename, const unsigned int I        
) {
#define daxpy(n, a, x, y) \
    for (int k = 0; k < (n); k++) {\
        (y)[k] += (a) * (x)[k]; \
    }

    if (dt == 0.) {
        fprintf(stderr, "Error: dt must be greater than 0\n");
        return -1;
    }

    if (gamma < 0. || gamma > 1.) {
        fprintf(stderr, "Error: gamma must be in [0, 1]\n");
        return -1;
    }

    if (beta < 0. || beta > 0.5) {
        fprintf(stderr, "Error: beta must be in [0, 0.5]\n");
        return -1;
    }

    if (Kbd->n != Mbd->n) {
        fprintf(stderr, "Error: K and M matrices must have the same size\n");
        return -1;
    }
    
    if (Kbd->k != Mbd->k) {
        fprintf(stderr, "Error: K and M matrices must have the same band size\n");
        return -1;
    }
    
    if (Kbd->n != 2*n) {
        fprintf(stderr, "Error: Dimension mismatch between the matrices and the state vectors\n");
        return -1;
    }

    if (I >= n) {
        fprintf(stderr, "Error: I must be a valid node index (I in [0, n-1])\n");
        return -1;
    }


    FILE *final = fopen(final_filename, "w");
    if (final == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", final_filename);
        return -1;
    }
    
    FILE *time = fopen(time_filename, "w");
    if (time == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", time_filename);
        fclose(final);
        return -1;
    }

    #if DEV_3_LOG == 1
    FILE *log = fopen("log.txt", "w");
    if (log == NULL) {
        fprintf(stderr, "Error: Could not open file log.txt\n");
        fclose(final);
        fclose(time);
        return -1;
    }
    #endif

    fprintf(time, "%.15le %.15le %.15le %.15le %.15le\n", 0., u[2*I], u[2*I+1], v[2*I], v[2*I+1]);

    int return_code = 0;

    CSRMatrix *M_csr = band_to_sym_csr(Mbd);
    CSRMatrix *K_csr = band_to_sym_csr(Kbd);

    SymBandMatrix *Mbd_cpy = allocate_sym_band_matrix(Mbd->n, Mbd->k);
    memcpy(Mbd_cpy->data, Mbd->data, sizeof(double) * Mbd->n * (Mbd->k + 1));
    sym_band_LDL(Mbd_cpy->data, Mbd_cpy->n, Mbd_cpy->k);

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
    CSRMatrix *rhs_csr = band_to_sym_csr(rhs_bd);

    double *p = (double *) malloc(2 * n * sizeof(double));
    double *rhs = (double *) malloc(2 * n * sizeof(double));
    double *interp = (double *) malloc(2 * n * sizeof(double));
    double *tmp = (double *) malloc(2 * n * sizeof(double));
    if (p == NULL || rhs == NULL || interp == NULL || tmp == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return_code = -1;
        goto free_matrices;
    }
    
    // p = Mv
    Matvec(M_csr->n, M_csr->row_ptr, M_csr->col_idx, M_csr->data, v, p);

    double t = 0.;
    while (t < T-dt/2.) {
        // Compute the right-hand side of first equation
        Matvec(rhs_csr->n, rhs_csr->row_ptr, rhs_csr->col_idx, rhs_csr->data, u, rhs);
        daxpy(2 * n, dt, p, rhs);
        
        // Solve the system, new u is in rhs 
        solve_sym_band(lhs_bd->data, lhs_bd->n, lhs_bd->k, rhs);
        
        // Linear interpolation of u
        memset(interp, 0, 2 * n * sizeof(double));
        daxpy(2 * n, (1-gamma),  u, interp);
        daxpy(2 * n,    gamma, rhs, interp);
        
        Matvec(K_csr->n, K_csr->row_ptr, K_csr->col_idx, K_csr->data, interp, tmp);
        
        // Update velocity
        daxpy(2 * n, -dt, tmp, p);
        
        memcpy(u, rhs, 2 * n * sizeof(double));
        
        memcpy(v, p, 2 * n * sizeof(double));
        solve_sym_band(Mbd_cpy->data, Mbd_cpy->n, Mbd_cpy->k, v);
        
        t += dt;

        fprintf(time, "%.15le %.15le %.15le %.15le %.15le\n", t, u[2*I], u[2*I+1], v[2*I], v[2*I+1]);
        
        #if DEV_3_LOG == 1

        Matvec(K_csr->n, K_csr->row_ptr, K_csr->col_idx, K_csr->data, u, tmp); 
        double Ep = 0.;
        for (int k = 0; k < 2 * n; k++)
            Ep += u[k] * tmp[k];
        Ep /= 2.;

        Matvec(M_csr->n, M_csr->row_ptr, M_csr->col_idx, M_csr->data, v, tmp);
        double Ek = 0.;
        for (int k = 0; k < 2 * n; k++)
            Ek += v[k] * tmp[k];
        Ek /= 2.;

        // The total energy should be constant
        fprintf(log, "%.15le %.15le %.15le %.15le\n", t, Ep, Ek, Ep + Ek);

        #endif
    }

    for (int k = 0; k < n; k++)
        fprintf(final, "%.15le %.15le %.15le %.15le\n", u[2*k], u[2*k+1], v[2*k], v[2*k+1]);

free_all:
    free(p);
    free(rhs);
    free(interp);
    free(tmp);
free_matrices:
    free_csr(M_csr);
    free_csr(K_csr);
    free_csr(rhs_csr);
    free_sym_band_matrix(lhs_bd);
    free_sym_band_matrix(rhs_bd);
    free_sym_band_matrix(Mbd_cpy);

    fclose(final);
    fclose(time);
    #if DEV_3_LOG == 1
    fclose(log);
    #endif

    return return_code;

#undef daxpy
}

#undef DEV_3_LOG