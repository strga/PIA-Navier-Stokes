#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

void writeDataToFile(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file '" << filename << "' for writing." << std::endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();

    if (file.good()) {
        std::cout << "Data successfully written to '" << filename << "'." << std::endl;
    } else {
        std::cerr << "Error occurred when writing to file '" << filename << "'." << std::endl;
    }
}

// Constants
const int N = 100;                       // Grid size
const double dom_length = 1.0;           // Domain length
const double dh = dom_length / (N - 1);  // Grid spacing
const double nu = 0.01;                  // Kinematic viscosity
const double URF = 0.8;                  // Under-relaxation factor for velocities
const double URF_p = 0.8;                // Under-relaxation factor for pressure

int main() {

    // Initialize variables
    std::vector<std::vector<double>> u(N + 1, std::vector<double>(N, 0.0)),
            v(N, std::vector<double>(N + 1, 0.0)),             //matrix with N rows and N+1 columns with zeros
            p(N + 1, std::vector<double>(N + 1, 1.0)),      //matrix with N+1 rows and N+1 columns with ones
            u_stag(N + 1, std::vector<double>(N, 0.0)),        //matrix with N+1 rows and N columns with zeros
            v_stag(N, std::vector<double>(N + 1, 0.0)),        //matrix with N rows and N+1 columns with zeros
            p_stag(N + 1, std::vector<double>(N + 1, 1.0)), //matrix with N rows and N+1 columns with zeros
            pc(N + 1, std::vector<double>(N + 1, 0.0)),     //matrix with N+1 rows and N+1 columns with zeros
            p_new(N + 1, std::vector<double>(N + 1, 1.0)),  //matrix with N+1 rows and N+1 columns with zeros
            u_new(N + 1, std::vector<double>(N, 0.0)),         //matrix with N+1 rows and N+1 columns with zeros
            v_new(N, std::vector<double>(N + 1, 0.0)),         //matrix with N rows and N+1 columns with zeros
            d_e(N + 1, std::vector<double>(N, 0.0)),           //matrix with N+1 rows and N columns with zeros
            d_n(N, std::vector<double>(N + 1, 0.0)),           //matrix with N rows and N+1 columns with zeros
            b(N + 1, std::vector<double>(N + 1, 0.0)),      //matrix with N+1 rows and N+1 columns with zeros
            u_final(N, std::vector<double>(N, 0.0)),              //Collocated final velocities
            v_final(N, std::vector<double>(N, 0.0)),              //Collocated final velocities
            p_final(N, std::vector<double>(N, 1.0));              //Collocated final pressure


    double err = 1.0; // Error for convergence check
    int iter = 0; // Iteration counter

    //for (int j = 0; j < N; ++j) {
    //    u_final[0][j] = 1.0;
    //}

    // Lid velocity initialization
    for (int j = 0; j < N; ++j) {
        u[0][j] = 2.0; // Top lid moving with double the unit velocity adjusted from 2.0
    }

    //err = 1.0; //adjusted from 0.0

    do {
        // Reset error for each iteration
        //err = 1.0; //adjusted from 0.0

        // Resetting the pressure-correction term pc
        //for (int i = 0; i <= N; ++i) {
        //    std::fill(pc[i].begin(), pc[i].end(), 0.0);
        //}


        // X-momentum equations - For Interior Control Volumes (CVs)
        // X-mom equations - For Interior CVs
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

        // X-momentum equations - For Boundary CVs
        for (int j = 0; j < N; ++j) {
            u_stag[0][j] = 2.0 - u_stag[1][j];
            u_stag[N][j] = -u_stag[N - 1][j]; //Maybe adjust here
        }
        for (int i = 1; i < N; ++i) {
            u_stag[i][0] = 0;
            u_stag[i][N - 1] = 0;
        }


        // Y-momentum equations - For Interior CVs
        // Y-momentum equations - For Interior CVs
        for (int i = 1; i < N-1; ++i) {
            for (int j = 1; j < N; ++j) { //<=
                double u_e = 0.5 * (u[i][j] + u[i + 1][j]);
                double u_w = 0.5 * (u[i][j - 1] + u[i + 1][j - 1]);
                double v_n = 0.5 * (v[i-1][j] + v[i][j]);
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


        // Boundary Conditions for u and v

        // Y-momentum equations - For Boundary CVs
        for (int i = 0; i < N; ++i) {
            v_stag[i][0] = -v_stag[i][1];
            v_stag[i][N] = -v_stag[i][N-1];
        }
        for (int j = 1; j < N - 1; ++j) { //
            v_stag[0][j] = 0; //v_stag[0][j] = 0;
            v_stag[N-1][j] = 0; //v_stag[N - 1][j] = 0;
        }

        // Pressure-correction term
        for (int i = 0; i <= N; ++i) {
            std::fill(pc[i].begin(), pc[i].end(), 0.0);
        }

        // Resetting the pressure-correction term pc
        for (int i = 0; i <= N; ++i) {
            std::fill(pc[i].begin(), pc[i].end(), 0.0);
        }

        // Solving for Continuity equation to get Pressure Correction - For Interior CVs
        // Solving for Continuity equation to get P - For Interior CVs
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                double a_E = -d_e[i][j] * dh;
                double a_W = -d_e[i][j - 1] * dh;
                double a_N = -d_n[i - 1][j] * dh;
                double a_S = -d_n[i][j] * dh;
                double a_P = a_E + a_W + a_N + a_S;
                b[i][j] = -(u_stag[i][j] - u_stag[i][j - 1]) * dh + (v_stag[i][j] - v_stag[i - 1][j]) * dh;

                pc[i][j] = (a_E * pc[i][j + 1] + a_W * pc[i][j - 1] + a_N * pc[i - 1][j] + a_S * pc[i + 1][j] + b[i][j]) / a_P;
            }
        }


        // Correcting the Pressure field and velocities
        // Correcting the Pressure field and velocities
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                p_new[i][j] = p[i][j] + URF_p * pc[i][j];
                //u_new[i][j] = u_stag[i][j] + URF * d_e[i][j] * (pc[i][j + 1] - pc[i][j]);
                //v_new[i][j] = v_stag[i][j] + URF * d_n[i][j] * (pc[i][j] - pc[i + 1][j]);
            }
        }


        // Continuity equation - Boundary Conditions for Pressure
        // Continuity equation - Boundary CVs for pressure
        for (int j = 0; j <= N; ++j) {
            p_new[0][j] = p_new[1][j];
            p_new[N][j] = p_new[N - 1][j];
        }
        for (int i = 0; i <= N; ++i) {
            p_new[i][0] = p_new[i][1];
            p_new[i][N] = p_new[i][N - 1];
        }
//+++++++++++++++++++++++++++++++++++++++++++++++ADD_HERE++++++++++++++++++++++++++++++++++

        // Correcting the velocities using pressure correction
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u_new[i][j] = u_stag[i][j] + URF * d_e[i][j] * (pc[i][j + 1] - pc[i][j]);
            }
        }

        // Applying boundary conditions for the corrected velocities
        for (int j = 0; j < N; ++j) {
            u_new[0][j] = 2.0 - u_new[1][j];
            u_new[N][j] = -u_new[N - 1][j];
        }
        for (int i = 1; i < N - 1; ++i) {
            u_new[i][0] = 0;
            u_new[i][N - 1] = 0;
        }

        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N; ++j) {
                v_new[i][j] = v_stag[i][j] + URF * d_n[i][j] * (pc[i][j] - pc[i + 1][j]);
            }
        }

        // Applying boundary conditions for the corrected velocities

        for (int i = 0; i < N; ++i) {
            v_new[i][0] = -v_new[i][1];
            v_new[i][N] = -v_new[i][N - 1];
        }
        for (int j = 1; j < N - 1; ++j) {
            v_new[0][j] = 0;
            v_new[N - 1][j] = 0;
        }


//+++++++++++++++++++++++++++++++++++++++++++++++ADD_HERE++++++++++++++++++++++++++++++++++
       // for (int j = 0; j < N; ++j) {
           // u_new[0][j] = 1.0;  // Set top lid velocity to 1.0 after correction
       // }

        // Update error value from residual term 'b' // i and j from for adjusted from i,j = 1 to i,j = 2
        err = 0;
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                err += std::abs(b[i][j]);
            }
        }

        // Update the values of variables for the next iteration
        u = u_new;
        v = v_new;
        p = p_new;
        //std::cout << "Iteration: " << iter << " Error: " << err << std::endl;
        ++iter;

    } while (err > 1e-6); // Convergence check


    // Populate u_final and v_final with computed data
    // Mapping staggered variables to collocated variables
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u_final[i][j] = 0.5 * (u[i][j] + u[i + 1][j]);
            v_final[i][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            p_final[i][j] = p[i][j]; // No change for pressure as its values are stored at a collocated grid
        }
    }


    // Check if the vectors are populated
    if (!u_final.empty() && !u_final[0].empty()) {
        writeDataToFile(u_final, "u_data.csv");
    } else {
        std::cout << "u_final is empty, skipping file write." << std::endl;
    }

    if (!v_final.empty() && !v_final[0].empty()) {
        writeDataToFile(v_final, "v_data.csv");
    } else {
        std::cout << "v_final is empty, skipping file write." << std::endl;
    }

    return 0;
}
