#pragma once

#include "model.h"
#include "utils.h"
#include "devoir_2.h"
#include "gmshc.h"
#include "utils_gmsh.h"

#include <string.h>
#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <time.h>


/**
 * @brief Reads the initial conditions from a file into a position and velocity arrays.
 * @param filename The name of the file containing the initial conditions.
 * @param u Pointer to the array for storing the position values.
 * @param v Pointer to the array for storing the velocity values.
 * @return The number of lines read from the file, or -1 on error.
 * @note The file should contain lines in the format: <ux_i> <uy_i> <vx_i> <vy_i>.
 *       The function allocates memory for the arrays and the caller is responsible for freeing it.
 *       The arrays will be of size 2*n, where n is the number of lines in the file.
 */
int read_initial_conditions(const char *filename, double **u, double **v);

/**
 * @brief Performs the Newmark algorithm on the given system.
 * @param Kbd Pointer to the stiffness matrix.
 * @param Mbd Pointer to the mass matrix.
 * @param u Pointer to the array for storing the position values.
 * @param v Pointer to the array for storing the velocity values.
 * @param n The number of nodes in the system.
 * @param dt The time step for the simulation.
 * @param T The total time for the simulation.
 * @param gamma The Newmark parameter for the velocity update.
 * @param beta The Newmark parameter for the position update.
 * @param final_filename The name of the file to save the final positions and velocities.
 * @param time_filename The name of the file to save position and velocity of a node I at each time step.
 * @param I The index of the node to save in time_filename.
 * @return 0 on success, or -1 on error.
 * @note The final position and velocity values will be stored in the arrays u and v.
 */
int newmark(
    FE_Model *model, double *u, double *v,
    const unsigned int n, const double dt, const double T, const double gamma, const double beta,
    const char *final_filename, const char *time_filename, const unsigned int I      
);
