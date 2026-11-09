#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <omp.h>

#include "Simulation.h"
#include "write_to_file.h"

using namespace std;

void runSimulation(int N, double domainLength, double kinematicViscosity)
{
    // Constants for the simulation
    const double dh = domainLength / (N - 1); // Grid spacing
    const double URF = 0.7;                  // nder-relaxation factor for velocities
    const double URF_p = 0.3;                // Under-relaxation factor for pressure
    const double nu = kinematicViscosity;    // Kinematic viscosity

    // Initialize matrices for simulation variables
    vector<vector<double>> u(N + 1, vector<double>(N, 0.0)), // Velocity in x-direction
    v(N, vector<double>(N + 1, 0.0)), // Velocity in y-direction
    p(N + 1, vector<double>(N + 1, 1.0)), // Pressure
    u_stag(N + 1, vector<double>(N, 0.0)), // Staggered grid for u
    v_stag(N, vector<double>(N + 1, 0.0)), // Staggered grid for v
    pc(N + 1, vector<double>(N + 1, 0.0)), // Pressure correction
    p_new(N + 1, vector<double>(N + 1, 1.0)), // New pressure values
    u_new(N + 1, vector<double>(N, 0.0)), // New u values
    v_new(N, vector<double>(N + 1, 0.0)), // New v values
    d_e(N + 1, vector<double>(N, 0.0)), // Coefficient for discretization
    d_n(N, vector<double>(N + 1, 0.0)), // Coefficient for discretization
    b(N + 1, vector<double>(N + 1, 0.0)), // Source term for pressure correction
    u_final(N, vector<double>(N, 0.0)), // Final u velocity
    v_final(N, vector<double>(N, 0.0)), // Final v velocity
    p_final(N, vector<double>(N, 1.0)); // Final pressure

    double err = 1.0; // Error for convergence
    int iter = 0;     // Iteration counter

    // Initialize top lid moving velocity (in both buffers, since u and u_new are swapped each iteration)
    for (int j = 0; j < N; ++j) {
        u_final[0][j] = 1.0;
        u_new[0][j] = 1.0;
        u[0][j] = 1.0;
    }

    do {
        // X-momentum: each (i,j) writes only its own cell and reads old u, v, p -> no races
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                double u_e = 0.5 * (u[i][j] + u[i][j + 1]);
                double u_w = 0.5 * (u[i][j] + u[i][j - 1]);
                double v_n = 0.5 * (v[i - 1][j] + v[i - 1][j + 1]);
                double v_s = 0.5 * (v[i][j] + v[i][j + 1]);

                double a_E = -0.5 * u_e * dh + nu;
                double a_W = 0.5 * u_w * dh + nu;
                double a_N = -0.5 * v_n * dh + nu;
                double a_S = 0.5 * v_s * dh + nu;
                double a_e = 0.5 * u_e * dh - 0.5 * u_w * dh + 0.5 * v_n * dh - 0.5 * v_s * dh + 4 * nu;

                double A_e = -dh;
                d_e[i][j] = A_e / a_e;

                u_stag[i][j] = (a_E * u[i][j + 1] + a_W * u[i][j - 1] + a_N * u[i - 1][j] + a_S * u[i + 1][j]) / a_e + d_e[i][j] * (p[i][j + 1] - p[i][j]);
            }
        }

        // Boundary conditions for u_stag (O(N), left serial)
        for (int j = 0; j < N; ++j) {
            u_stag[0][j] = 1.0 - u_stag[1][j]; //top
            u_stag[N][j] = -u_stag[N - 1][j]; //bottom
        }
        for (int i = 1; i < N; ++i) {
            u_stag[i][0] = 0; //left
            u_stag[i][N - 1] = 0; //right
        }

        // Y-momentum
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N; ++j) {
                double u_e = 0.5 * (u[i][j] + u[i + 1][j]);
                double u_w = 0.5 * (u[i][j - 1] + u[i + 1][j - 1]);
                double v_n = 0.5 * (v[i - 1][j] + v[i][j]);
                double v_s = 0.5 * (v[i][j] + v[i + 1][j]);

                double a_E = -0.5 * u_e * dh + nu;
                double a_W = 0.5 * u_w * dh + nu;
                double a_N = -0.5 * v_n * dh + nu;
                double a_S = 0.5 * v_s * dh + nu;
                double a_n = 0.5 * u_e * dh - 0.5 * u_w * dh + 0.5 * v_n * dh - 0.5 * v_s * dh + 4 * nu;

                double A_n = -dh;
                d_n[i][j] = A_n / a_n;

                v_stag[i][j] = (a_E * v[i][j + 1] + a_W * v[i][j - 1] + a_N * v[i - 1][j] + a_S * v[i + 1][j]) / a_n + d_n[i][j] * (p[i][j] - p[i + 1][j]);
            }
        }

        // Boundary conditions for v_stag
        for (int i = 0; i < N; ++i) {
            v_stag[i][0] = -v_stag[i][1]; //left
            v_stag[i][N] = -v_stag[i][N - 1]; //right
        }
        for (int j = 1; j < N - 1; ++j) {
            v_stag[0][j] = 0; //top
            v_stag[N - 1][j] = 0; //bottom
        }

        // Reset pressure correction
#pragma omp parallel for
        for (int i = 0; i <= N; ++i) {
            std::fill(pc[i].begin(), pc[i].end(), 0.0);
        }

        // Continuity imbalance (source term), computed once and in parallel
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                b[i][j] = -(u_stag[i][j] - u_stag[i][j - 1]) * dh + (v_stag[i][j] - v_stag[i - 1][j]) * dh;
            }
        }

        // Pressure correction: red-black Gauss-Seidel sweep
        for (int colour = 0; colour < 2; ++colour) {
#pragma omp parallel for
            for (int i = 1; i < N; ++i) {
                int jStart = 1 + ((i + 1 + colour) % 2);
                for (int j = jStart; j < N; j += 2) {
                    double a_E = -d_e[i][j] * dh;
                    double a_W = -d_e[i][j - 1] * dh;
                    double a_N = -d_n[i - 1][j] * dh;
                    double a_S = -d_n[i][j] * dh;
                    double a_P = a_E + a_W + a_N + a_S;

                    pc[i][j] = (a_E * pc[i][j + 1] + a_W * pc[i][j - 1] + a_N * pc[i - 1][j] + a_S * pc[i + 1][j] + b[i][j]) / a_P;
                }
            }
        }

        // Correct the pressure field
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                p_new[i][j] = p[i][j] + URF_p * pc[i][j];
            }
        }

        // Boundary conditions for pressure
        for (int j = 0; j <= N; ++j) {
            p_new[0][j] = p_new[1][j];
            p_new[N][j] = p_new[N - 1][j];
        }
        for (int i = 0; i <= N; ++i) {
            p_new[i][0] = p_new[i][1];
            p_new[i][N] = p_new[i][N - 1];
        }

        // Correct u with the pressure correction
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u_new[i][j] = u_stag[i][j] + URF * d_e[i][j] * (pc[i][j + 1] - pc[i][j]);
            }
        }

        // Boundary conditions for corrected u
        for (int j = 0; j < N; ++j) {
            u_new[0][j] = 1.0 - u_new[1][j];
            u_new[N][j] = -u_new[N - 1][j];
        }
        for (int i = 1; i < N - 1; ++i) {
            u_new[i][0] = 0;
            u_new[i][N - 1] = 0;
        }

        // Correct v with the pressure correction
#pragma omp parallel for collapse(2)
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N; ++j) {
                v_new[i][j] = v_stag[i][j] + URF * d_n[i][j] * (pc[i][j] - pc[i + 1][j]);
            }
        }

        // Boundary conditions for corrected v
        for (int i = 0; i < N; ++i) {
            v_new[i][0] = -v_new[i][1];
            v_new[i][N] = -v_new[i][N - 1];
        }
        for (int j = 1; j < N - 1; ++j) {
            v_new[0][j] = 0;
            v_new[N - 1][j] = 0;
        }

        // Residual from the continuity imbalance
        err = 0;
#pragma omp parallel for collapse(2) reduction(+:err)
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                err += std::abs(b[i][j]);
            }
        }

        std::swap(u, u_new);
        std::swap(v, v_new);
        std::swap(p, p_new);
        ++iter;

    } while (err > 1e-6);

    cout << "Iterations: " << iter << endl;

    // Mapping staggered variables to collocated variables
#pragma omp parallel for collapse(2)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u_final[i][j] = 0.5 * (u[i][j] + u[i + 1][j]);
            v_final[i][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            p_final[i][j] = p[i][j];
        }
    }

    // Write final data to .csv
    if (!u_final.empty() && !u_final[0].empty()) {
        writeDataToFile(u_final, "u_data.csv");
    } else {
        cout << "u_final is empty, skipping file write." << endl;
    }

    if (!v_final.empty() && !v_final[0].empty()) {
        writeDataToFile(v_final, "v_data.csv");
    } else {
        cout << "v_final is empty, skipping file write." << endl;
    }
}
