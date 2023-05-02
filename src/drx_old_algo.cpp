#include <iostream>
#include <random>
#include <chrono>
#include "helpers.cpp"
#include "drx_grid.cpp"

int main() {
    // Invoking GnuPlot in the background.
    system("gnuplot cadrx.gp &");

    // seeding random number initialization with current time
    srand((unsigned)time(0));

    _grid_ grid;
    _grid_ u_grid;
    deep_copy_grid(&grid, &u_grid);
    // print_grid(grid.cell);
    // cout << "\n\n";
    // print_array(grid.grain_num);

    // timing variables
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();;

    // bool on_border = 0;
    // int x = 0, y = 0;
    // int nx = 0, ny = 0;
    // int x_start = 0, x_end = 3, y_start = 0, y_end = 3;

    // float cell_gamma = 0;
    // float delta_p = 0;
    // float misorientation = 0;
    // float grain_size = 0;
    
    float eps = 0;
    float delta_t = 0;
    float delta_eps = 0;
    float delta_n = 0;

    // vector<int> potential_nucleus_x {-1};
    // vector<int> potential_nucleus_y {-1};
    // int pot_n[GRID_SIZE][GRID_SIZE];

    float p_nucleation = 0;
    float p_growth = 0;
    float p_random = 0;
    // int pn_i = 0;
    // int nucleus_counter = 0;
    int neighbour_info = 0;

    int iteration = 0;

    int i = 0, j = 0;

    write_to_file(grid.cell);

    while (eps < EPS_FINAL) {
        cout << "Iteration: " << ++iteration << endl;
        delta_t = CELL_SIZE / grid.v_max;
        delta_eps = EPS_RATE * delta_t;
        eps += delta_eps;
        if (eps > EPS_CR) {
            ff (i, 0, GRID_SIZE) {
                ff(j, 0, GRID_SIZE) {
                    if (grid.cell[i][j].N_recrystallized == 1) {
                        if(!grid.propagate_grain_boundary(i, j, iteration, &u_grid)) {
                            u_grid.cell[i][j].update_dislocation_density(delta_eps);
                        }
                    } else if (grid.potential_nucleus(i, j)) {
                        p_random = (float)rand() / (float)RAND_MAX;
                        p_nucleation = 0.025 * grid.cell[i][j].dislocation_density / grid.p_max;
                        if (p_random < p_nucleation) {
                            u_grid.cell[i][j].nucleate();
                        } else {
                            u_grid.cell[i][j].update_dislocation_density(delta_eps);
                        }
                    } else {
                        u_grid.cell[i][j].update_dislocation_density(delta_eps);
                    }
                }
            }
            write_to_file(u_grid.cell);
            cin.get();
        } else {
            ff (i, 0, GRID_SIZE) {
                ff (j, 0, GRID_SIZE) {
                    u_grid.cell[i][j].update_dislocation_density(delta_eps);
                }
            }
        }
        deep_copy_grid(&u_grid, &grid);

        grid.average_p();
        cout << "Average p: " << grid.p_avg << endl;
        
    }
    write_to_file(grid.cell);

    return 0;
}