#include <iostream>
#include <random>
#include <chrono>
#include "helpers.cpp"
#include "drx_grid.cpp"

int main() {
    // Invoking GnuPlot in the background.
    // system("gnuplot cadrx.gp &");

    // seeding random number initialization with current time
    srand((unsigned)time(0));

    _grid_ grid;
    _grid_ u_grid;
    deep_copy_grid(&grid, &u_grid);
    
    float eps = 0;
    float delta_t = 0;
    float delta_eps = 0;
    float delta_n = 0;

    vector<int> potential_nucleus_x {-1};
    vector<int> potential_nucleus_y {-1};
    int nucleated[GRID_SIZE][GRID_SIZE];
    float p_nucleation = 0;
    float p_growth = 0;
    float p_random = 0;
    int pn_i = 0;
    int nucleus_counter = 0;

    int iteration = 0;

    float xdrx = 0;
    float stress = 0;

    ofstream stress_strain;
    stress_strain.open("strainVstress.dat");

    write_to_file(grid.cell);

    while (eps < EPS_FINAL) {
        // cout << "Iteration: " << ++iteration << endl;
        delta_t = CELL_SIZE / grid.v_max;
        delta_eps = EPS_RATE * delta_t;
        eps += delta_eps;
        if (eps > EPS_CR) {
            pn_i = 0;
            nucleus_counter = 0;

            // Calculate number of new nucleus to be formed
            delta_n = (int)(nucleation_rate * delta_t);

            ff(i, 0, GRID_SIZE) {
                ff(j, 0, GRID_SIZE) {
                    nucleated[i][j] = 0;
                    if (grid.potential_nucleus(i, j)) {
                        potential_nucleus_x.push_back(i);
                        potential_nucleus_y.push_back(j);
                        nucleus_counter++;
                    }
                }
            }

            if (nucleus_counter >= delta_n) {
                ff(i, 0, delta_n) {
                    pn_i = rand() % potential_nucleus_x.size();
                    if (potential_nucleus_y[pn_i] == -1) {
                        i--;
                        continue;
                    }
                    if (u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized == 1) {
                        i--;
                        continue;
                    }
                    p_random = rand() / RAND_MAX;
                    p_nucleation = grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].dislocation_density / grid.p_max;
                    if (p_random < p_nucleation && grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized != 1) {
                        u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].nucleate();
                        grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].nucleate();
                        potential_nucleus_x[pn_i] = potential_nucleus_y[pn_i] = -1;
                        nucleated[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]] = 1;
                    } 
                    // else {
                    //     u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].update_dislocation_density(delta_eps);
                    // }
                }
            }

            nucleus_counter = 0;
            
            ff (i, 0, GRID_SIZE) {
                ff(j, 0, GRID_SIZE) {
                    if (grid.cell[i][j].N_recrystallized == 1) {
                        nucleus_counter++;
                        // grid.propagate_grain_boundary(i, j, iteration, &u_grid);
                        // u_grid.cell[i][j].update_dislocation_density(delta_eps);
                        if(!grid.propagate_grain_boundary(i, j, iteration, &u_grid)) {
                            // u_grid.cell[i][j].update_dislocation_density(delta_eps);
                        }
                    } else if (nucleated[i][j] != 1) {
                        u_grid.cell[i][j].update_dislocation_density(delta_eps);
                    } else {
                        u_grid.cell[i][j].update_dislocation_density(delta_eps);
                    }
                }
            }
            // write_to_file(u_grid.cell);
            // cin.get();
        } else {
            ff (i, 0, GRID_SIZE) {
                ff (j, 0, GRID_SIZE) {
                    u_grid.cell[i][j].update_dislocation_density(delta_eps);
                }
            }
        }
        deep_copy_grid(&u_grid, &grid);
        /* DO NOT DO ANYTHING ABOVE THIS */

        xdrx = (float)nucleus_counter / (GRID_SIZE * GRID_SIZE);
        // cout << "xDrx: " << xdrx << endl;
        grid.average_p();
        // cout << "P Average: " << grid.p_avg << endl;
        stress = alpha * mu * b * sqrt(grid.p_avg) / 1000000;
        // cout << "Mu: " << mu << endl;
        stress_strain << eps << "\t" << (stress) << "\n";
        // cout << "Stress: " << stress << endl;
        potential_nucleus_x.clear();
        potential_nucleus_y.clear();
        // write_to_file(grid.cell);

        // cout << "Nucleus Counter: " << nucleus_counter << endl;

        // stop iterations if all the cells have nucleated
        if (nucleus_counter == (GRID_SIZE * GRID_SIZE)) {
            break;
        }
    }
    stress_strain.close();
    write_to_file(grid.cell);

    system("gnuplot ../data/strainVSstress.gp -p &");

    return 0;
}
