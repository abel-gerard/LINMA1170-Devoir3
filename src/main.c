#include "devoir_2.h"
#include "devoir_3.h"
#include "utils.h"
#include "model.h"
#include "utils_gmsh.h"
#include <math.h>
#include <cblas.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERBOSE 1
#define PRECISION 10

void display_info(FE_Model *model, int step, struct timespec ts[4]) {

    char *m_str[3] = {"Plane stress", "Plane strain", "Axisymmetric"};
    char *r_str[4] = {"No", "X", "Y", "RCMK"};

    if (step == 1) {
        printf(
            "\n===========  Linear elasticity simulation - FEM  ===========\n\n"
        );
        printf("%30s = %s\n", "Model", model->model_name);
        printf("%30s = %s\n", "Model type", m_str[model->m_type]);
        printf("%30s = %.3e\n", "Young's Modulus E", model->E);
        printf("%30s = %.3e\n", "Poisson ratio nu", model->nu);
        printf("%30s = %.3e\n", "Density rho", model->rho);
        printf("%30s = %.3e\n\n", "Frequency", (1.875*1.875/(4.*M_PI*model->L_ref*model->L_ref))*sqrt(model->E*0.0956*0.0956/model->rho));
    } else if (step == 2) {
        char *e_str = (model->e_type == TRI) ? "Triangle" : "Quadrilateral";
        printf("%30s = %s\n", "Element type", e_str);
        printf("%30s = %zu\n", "Number of elements", model->n_elem);
        printf("%30s = %zu\n", "Number of nodes", model->n_node);
        printf("%30s = %s\n", "Renumbering", r_str[model->renum]);
        printf("%30s = %zu\n\n", "Matrix bandwidth", 2 * model->node_band + 1);
    }
}

int main(int argc, char *argv[]) {

    int ierr;
    double mesh_size_ratio;
    if ((argc < 9) || (sscanf(argv[2], "%lf", &mesh_size_ratio)) != 1) {
        printf("Usage: \n./deformation <model> <mesh_size_ratio> <T> <dt> <initial.txt> <final.txt> <time.txt> <I>\n");
        printf("model: one of the model implemented in models/\n");
        printf("mesh_size_ratio: mesh size factor\n");
        return -1;
    }

    // For example "./deformation fork 1 1 0.01 initial_fork_1.0.txt final_fork.txt time_fork.txt 50"
    const double T = atof(argv[3]);
    const double dt = atof(argv[4]);
    const char *initial_file = argv[5];
    const char *final_file = argv[6];
    const char *time_file = argv[7];
    const int _I = atoi(argv[8]); 

    const double gamma = 0.5;
    const double beta = 0.25;

    // Simulation parameters
    const ElementType e_type = TRI;
    const Renumbering renum = RENUM_NO;  // let gmsh do the RCMK renumbering

    FE_Model *model = create_FE_Model(argv[1], e_type, renum);
    display_info(model, 1, NULL);

    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 2, &ierr);
    model->mesh_model(mesh_size_ratio, e_type);

    load_mesh(model);
    renumber_nodes(model);
    display_info(model, 2, NULL);
    assemble_system(model);
    double *rhs = (double *)calloc(2 * model->n_node, sizeof(double));
    double *sol = (double *)calloc(2 * model->n_node, sizeof(double));
    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);
    
    SymBandMatrix *Kbd = model->K;
    SymBandMatrix *Mbd = model->M;
    
    double *u;
    double *v;
    int n = read_initial_conditions(initial_file, &u, &v);
    if (n < 0) {
        fprintf(stderr, "Error: Could not read initial conditions\n");
        return -1;
    }

    if (newmark(
        model, u, v,
        n, dt, T, gamma, beta,
        final_file, time_file, _I
    ) < 0) {
        fprintf(stderr, "Error: Newmark algorithm failed\n");
        return -1;
    }

    // CSRMatrix *Ksp = band_to_sym_csr(Kbd);
    // CSRMatrix *Msp = band_to_sym_csr(Mbd);
    // double eps = 1e-8;
    // CG(Ksp->n, Ksp->row_ptr, Ksp->col_idx, Ksp->data, rhs, sol, eps);    
    // display_sol(model, sol);

    // display_sol(model, u, NULL);    

    // Free stuff
    free(u);
    free(v);

    // free_csr(Ksp);
    // free_csr(Msp);
    gmshFinalize(&ierr);
    free(sol);
    free(rhs);
    free_FE_Model(model);
    return 0;
}
