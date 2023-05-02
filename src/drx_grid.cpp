#ifndef _DRX_GRID_CPP_
#define _DRX_GRID_CPP_

#include "drx_grid.h"
#include <fstream>

void grain_cell::nucleate() {
    this -> dislocation_density = DIS_DEN_RESET;
    this -> N_recrystallized = 1;
    this -> grain_number = rand() % STATES + 1;
    this -> orientation = ((float)rand() / (float)RAND_MAX) * (float)MAX_ORIENTATION;
}

void grain_cell::update_dislocation_density(float delta_str) {
    float delta_p = (k1 * sqrt(this -> dislocation_density) - k2 * (this -> dislocation_density)) * delta_str;
    // if (this -> dislocation_density + delta_p <= DIS_DEN_MAX) {
    //     this -> dislocation_density += delta_p;
    // }
    this -> dislocation_density += delta_p;
}

int _grid_::check_neighbours(int i, int j) {
    int x = 0, y = 0;
    int nx = 0, ny = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;
    if (i == 0) x_start++;
    if (j == 0) y_start++;
    if (i == GRID_SIZE - 1) x_end--;
    if (j == GRID_SIZE - 1) y_end--;
    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y < y_end; y++) {
            nx = i + x - 1;
            ny = j + y - 1;
            if (i == x && j == y) continue;
            if (this -> cell[i][j].grain_number != this -> cell[nx][ny].grain_number) {
                return nx * pow(10, encoder + 1) + ny * pow(10, 1) + 1;
            }
        }
    }
    return 0;
}

void _grid_::voronoiInit(grain_cell mat[GRID_SIZE][GRID_SIZE]) {
    
}

void _grid_::monteCarloInit(grain_cell mat[GRID_SIZE][GRID_SIZE]) {
    int i, j;
    // int* arr = (int*)malloc(sizeof(int) * STATES);
    int arr[STATES];
    int index = 0, temp_max = 0;
    int not_similar_i = 0, not_similar_f = 0, delta = 0;
    float prob, r_no;
    for (i = 0; i < ITERATIONS_MC; i++) {
        for (j = 0; j < STATES; j++) {
            arr[j] = 0;
        }
        not_similar_i = 0;
        not_similar_f = 0;
        delta = 0;
        int x = rand() % GRID_SIZE;
        int y = rand() % GRID_SIZE;
        if ( x > 0 && y > 0) {
            arr[(int)mat[x - 1][y - 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x - 1][y - 1].grain_number) not_similar_i++;
        }
        if (x > 0) {
            arr[(int)mat[x - 1][y].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x - 1][y].grain_number) not_similar_i++;
        }
        if(x > 0 && y < GRID_SIZE - 1) {
            arr[(int)mat[x - 1][y + 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x - 1][y].grain_number) not_similar_i++;
        }
        if(y > 0) {
            arr[(int)mat[x][y - 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x][y - 1].grain_number) not_similar_i++;
        }
        if(y < GRID_SIZE - 1) {
            arr[(int)mat[x][y + 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x][y + 1].grain_number) not_similar_i++;
        }
        if(x < GRID_SIZE - 1 && y > 0) {
            arr[(int)mat[x + 1][y - 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x + 1][y - 1].grain_number) not_similar_i++;
        }
        if(x < GRID_SIZE - 1) {
            arr[(int)mat[x + 1][y].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x + 1][y - 1].grain_number) not_similar_i++;
        }
        if(x < GRID_SIZE - 1 && y < GRID_SIZE - 1) {
            arr[(int)mat[x + 1][y + 1].grain_number - 1] += 1;
            if (mat[x][y].grain_number != mat[x + 1][y + 1].grain_number) not_similar_i++;
        }
        temp_max = 0;
        for (j = 0; j < STATES; j++) {
            if (arr[j] > temp_max) {
                temp_max = arr[j];
                index = j + 1;
            }
        }
        not_similar_f = 8 - temp_max;
        delta = not_similar_f - not_similar_i;
        if (delta < 0) {
            mat[x][y].grain_number = index;
        } else {
            prob = abs(delta) / 6;
            r_no = (float) rand() / RAND_MAX;
            if (r_no < prob) {
                mat[x][y].grain_number = index;
            }
        }
    }

    ff(i, 0, 4) {
        int x_ = -1, y_ = -1;
        int x = 0, y = 0;
        if (i % 4 == 0) {
            x_ = 0;
        } else if (i % 4 == 1) {
            x_ = GRID_SIZE - 1;
        } else if (i % 4 == 2) {
            y_ = 0;
        } else if (i % 4 == 3) {
            y_ = GRID_SIZE - 1;
        }
        ff(k, 0, GRID_SIZE) {
            if (x_ == -1) {
                x = k;
                y = y_;
            }
            if (y_ == -1) {
                y = k;
                x = x_;
            }
            // mat[x][y].grain_number = 0;
            for (j = 0; j < STATES; j++) {
                arr[j] = 0;
            }
            not_similar_i = 0;
            not_similar_f = 0;
            delta = 0;
            if ( x > 0 && y > 0) {
                arr[(int)mat[x - 1][y - 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x - 1][y - 1].grain_number) not_similar_i++;
            }
            if (x > 0) {
                arr[(int)mat[x - 1][y].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x - 1][y].grain_number) not_similar_i++;
            }
            if(x > 0 && y < GRID_SIZE - 1) {
                arr[(int)mat[x - 1][y + 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x - 1][y].grain_number) not_similar_i++;
            }
            if(y > 0) {
                arr[(int)mat[x][y - 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x][y - 1].grain_number) not_similar_i++;
            }
            if(y < GRID_SIZE - 1) {
                arr[(int)mat[x][y + 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x][y + 1].grain_number) not_similar_i++;
            }
            if(x < GRID_SIZE - 1 && y > 0) {
                arr[(int)mat[x + 1][y - 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x + 1][y - 1].grain_number) not_similar_i++;
            }
            if(x < GRID_SIZE - 1) {
                arr[(int)mat[x + 1][y].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x + 1][y - 1].grain_number) not_similar_i++;
            }
            if(x < GRID_SIZE - 1 && y < GRID_SIZE - 1) {
                arr[(int)mat[x + 1][y + 1].grain_number - 1] += 1;
                if (mat[x][y].grain_number != mat[x + 1][y + 1].grain_number) not_similar_i++;
            }
            temp_max = 0;
            for (j = 0; j < STATES; j++) {
                if (arr[j] > temp_max) {
                    temp_max = arr[j];
                    index = j + 1;
                }
            }
            mat[x][y].grain_number = index;
            // not_similar_f = 8 - temp_max;
            // delta = not_similar_f - not_similar_i;
            // if (delta < 0) {
            //     mat[x][y].grain_number = index;
            // } else {
            //     prob = abs(delta) / 6;
            //     r_no = (float) rand() / RAND_MAX;
            //     if (r_no < prob) {
            //         mat[x][y].grain_number = index;
            //     }
            // }
        }
    }

}

float _grid_::calculate_cell_velocity(float delta_p, float misorientation, float grain_size) {
    int nx = 0, ny = 0;
    float cell_gamma = 0;
    float P = 0;

    cell_gamma = gamma_m * sin(2 * misorientation * M_PI / 180) * (1 - r_gamma * log(sin(2 * misorientation * M_PI / 180)));
    
    M = mobility(misorientation);
    P = abs(tau * delta_p - 2 * cell_gamma / grain_size); // In Liu & Hallberg

    return M * P;
}

float _grid_::calculate_velocities() {
    int x = 0, y = 0;
    int nx = 0, ny = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;
    float cell_gamma = 0;
    float P = 0;
    float delta_p = 0;
    float misorientation = 0;
    bool on_border = 0;
    float grain_size = 0.00001;
    float v_max = 0;
    int neighbour_info = 0;

    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            neighbour_info = 0;
            
            neighbour_info = this -> check_neighbours(i, j);

            if (neighbour_info != 0) {
                // cout << neighbour_info << endl;
                nx = (int)(neighbour_info / pow(10, encoder + 1));
                neighbour_info = neighbour_info % (int)pow(10, encoder + 1);
                ny = (int)neighbour_info / 10;
                // cout << "For " << nx << " and " << ny << endl;
                misorientation = abs(this -> cell[i][j].orientation - this -> cell[nx][ny].orientation);
                cell_gamma = gamma_m * sin(2 * misorientation * M_PI / 180) * (1 - r_gamma * log(sin(2 * misorientation * M_PI / 180)));
                // cell_gamma = gamma_l(misorientation);
                delta_p = abs(this -> cell[i][j].dislocation_density - this -> cell[nx][ny].dislocation_density);
                grain_size = sqrt(this -> grain_size_nc[this -> grain_num[i][j]] * CELL_SIZE * CELL_SIZE / M_PI);
                M = mobility(misorientation);
                P = abs(tau * delta_p - 2 * cell_gamma / grain_size); // In Liu & Hallberg
                // P = tau * delta_p; // In Popova

                // p_ciritcal = pow((20 * cell_gamma * EPS_RATE) / (3 * b * sqrt(this -> cell[i][j].dislocation_density) * M * tau * tau), 1.0 / 3.0);
                // cout << "critical p: " << p_ciritcal << endl;

                this -> cell[i][j].gb_velocity = M * P;

                if (this -> cell[i][j].gb_velocity > v_max) {
                    // cout <<  this -> cell[i][j].gb_velocity;
                    v_max = this -> cell[i][j].gb_velocity;
                }
            }
        }
    }
    // cout << "Finished Calculating Velocities!!\n";
    return v_max;
}

void _grid_::reset_grain_variables() {
    this -> grain_size_nc.clear();
    grain_size_nc = {0};
    ff(i, 0, GRID_SIZE) {
        // this -> grain_size_nc[i % STATES] = 0;
        ff(j, 0, GRID_SIZE) {
            this -> grain_num[i][j] = 0;
        }
    }
}

void _grid_::saturate_grain(int i, int j, int parent_grain_count, int parent_grain_number, bool seed) {
    if (i < 0 || j < 0 || i >= GRID_SIZE || j >= GRID_SIZE || (this -> cell[i][j].grain_number != parent_grain_number)) {
        return;
    }
    if (!seed) {
        if (this -> grain_num[i][j] != 0) {
            return;
        }
    }
    // cout << i << " " << j << " " << parent_grain_count << " " << parent_grain_number << "\n";
    this -> grain_num[i][j] = parent_grain_count;
    this -> grain_size_nc[parent_grain_count]++;

    int nx = 0, ny = 0;
    int x = 0, y = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;
    x = y = 0;
    nx = ny = 0;
    x_start = y_start = 0;
    x_end = y_end = 3;
    if (i == 0) x_start++;
    if (j == 0) y_start++;
    if (i == GRID_SIZE - 1) x_end--;
    if (j == GRID_SIZE - 1) y_end--;
    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y < y_end; y++) {
            if (i == nx && j == ny) continue;
            nx = i + x - 1;
            ny = j + y - 1;
            saturate_grain(nx, ny, parent_grain_count, parent_grain_number, 0);
        }
    }
    return;
}

void _grid_::set_grain_numbers() {
    reset_grain_variables();
    int x = 0, y = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;
    int nx = 0, ny = 0;
    int grain_counter = 0;
    this -> grain_num[0][0] = ++grain_counter;
    bool found_seed = 0;
    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            if (i == 0 && j == 0) {
                this -> grain_size_nc.push_back(0);
                saturate_grain(i, j, grain_counter, this -> cell[i][j].grain_number, 1);
                continue;
            }
            if (this -> grain_num[i][j] != 0) {
                continue;
            } else {
                this -> grain_num[i][j] = ++grain_counter;
                this -> grain_size_nc.push_back(0);
                saturate_grain(i, j, grain_counter, this -> cell[i][j].grain_number, 1);
            }
        }
    }
    // cout << "Total Grains: " << grain_counter << "\n";
    this -> grain_size_nc.pop_back();
    
    ff(i, 0, this -> grain_size_nc.size()) {
        if (i == this -> grain_size_nc.size() - 1) {
            this -> grain_size_nc.pop_back();
            break;
        }
        this -> grain_size_nc[i] = this -> grain_size_nc[i + 1];
        // cout << this -> grain_size_nc[i] << "\n";
    }
    // cout << "Grain Count Vector Size: " << grain_size_nc.size() << "\n";
}

void _grid_::average_p() {
    float total_p = 0.0;

    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            if (this -> cell[i][j].dislocation_density > this -> p_max) {
                this -> p_max = this -> cell[i][j].dislocation_density;
            }
            total_p += this -> cell[i][j].dislocation_density;
        }
    }
    this -> p_avg = (float)total_p / (float)(GRID_SIZE * GRID_SIZE);
}

bool _grid_::potential_nucleus(int i, int j) {
    if (this -> cell[i][j].N_recrystallized == 1) {
        return 0;
    }
    int neighbour_info = 0;
    int nx = 0, ny = 0;
    float misorientation = 0;
    float cell_gamma = 0;
    neighbour_info = this -> check_neighbours(i, j);
    if (neighbour_info != 0) {
        nx = (int)(neighbour_info / pow(10, encoder + 1));
        neighbour_info = neighbour_info % (int)pow(10, encoder + 1);
        ny = (int)neighbour_info / 10;
        misorientation = abs(this -> cell[i][j].orientation - this -> cell[nx][ny].orientation);
        cell_gamma = gamma_m * sin(2 * misorientation * M_PI / 180) * (1 - r_gamma * log(sin(2 * misorientation * M_PI / 180)));
        M = mobility(misorientation);
        p_ciritcal = pow((20 * cell_gamma * EPS_RATE) / (3 * b * sqrt(this -> cell[i][j].dislocation_density) * M * tau * tau), 1.0 / 3.0);
        if (misorientation > CRITICAL_MISORIENTATION && this -> cell[i][j].dislocation_density > p_ciritcal) {
            return 1;
        }
    }
    return 0;
}

bool _grid_::consume_recrystallized_nuclei(int i, int j, int p_i, int p_j) {
    int nx = 0, ny = 0;
    int x = 0, y = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;

    int not_similar_i = 0;
    int not_similar_f = 0;
    int delta = 0;

    if (i == 0) x_start++;
    if (j == 0) y_start++;
    if (i == GRID_SIZE - 1) x_end--;
    if (j == GRID_SIZE - 1) y_end--;
    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y < y_end; y++) {
            nx = i + x - 1;
            ny = j + y - 1;
            if (i == nx && j == ny) continue;
            if (this -> cell[i][j].grain_number != this -> cell[nx][ny].grain_number) {
                not_similar_i++;
            }

            if (this -> cell[p_i][p_j].grain_number != this -> cell[nx][ny].grain_number) {
                not_similar_f++;
            }
        }
    }

    delta = not_similar_f - not_similar_i;
    // cout << "delta: " << delta << endl;
    if (delta < 0) {
        return 1;
    } else {
        return 0;
    }
    // if (not_similar_i == (x_end - x_start) * (y_end - y_start) - 1) {
    //     return 1;
    // } else {
    //     return 0;
    // }
}

bool _grid_::propagate_grain_boundary(int i, int j, int k, _grid_ *update_grid) {
    int nx = 0, ny = 0;
    int x = 0, y = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;

    float misorientation = 0;
    float delta_p = 0;
    float grain_size = 0;
    float gb_vel = 0;

    float p_random = 0;
    float p_growth = 0;

    bool propagated = 0;

    if (i == 0) x_start++;
    if (j == 0) y_start++;
    if (i == GRID_SIZE - 1) x_end--;
    if (j == GRID_SIZE - 1) y_end--;
    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y < y_end; y++) {
            nx = i + x - 1;
            ny = j + y - 1;
            if (i == nx && j == ny) continue;
            if ((x == x_start && y == y_start) || (x == x_end - 1 && y == y_end - 1)) {
                if (k % 2 == 1) continue;
            }
            if ((x == x_start && y == y_end - 1) || (x == x_end - 1 && y == y_start)) {
                if (k % 2 == 0) continue;
            }
            if (this -> cell[nx][ny].N_recrystallized == 1) {
                if (this -> cell[i][j].grain_number == this -> cell[nx][ny].grain_number) continue;
                if (consume_recrystallized_nuclei(nx, ny, i, j)) {
                    update_grid -> cell[nx][ny].grain_number = this -> cell[i][j].grain_number;
                    update_grid -> cell[nx][ny].dislocation_density = this -> cell[i][j].dislocation_density;
                    update_grid -> cell[nx][ny].N_recrystallized = 1;
                    update_grid -> cell[nx][ny].orientation = this -> cell[i][j].orientation;
                    propagated = 1;
                }

            } else {
            // if (this -> cell[nx][ny].N_recrystallized == 0) {
                misorientation = abs(this -> cell[i][j].orientation - this -> cell[nx][ny].orientation);
                delta_p = abs(this -> cell[i][j].dislocation_density - this -> cell[nx][ny].dislocation_density);
                grain_size = sqrt(this -> grain_size_nc[this -> grain_num[i][j]] * CELL_SIZE * CELL_SIZE / M_PI);
                gb_vel = this -> calculate_cell_velocity(delta_p, misorientation, grain_size);
                if (gb_vel > update_grid -> v_max) {
                    update_grid -> v_max = gb_vel;
                }
                p_growth = gb_vel / this -> v_max;
                p_random = (float)rand() / (float)RAND_MAX;
                // The function in the second contition is used here because it does not demand any of the cells
                // to be already nucleated. The name suggests that the nuclei to be consumed should be already nucleated but the implementation does not demand that.
                if (p_random < p_growth || consume_recrystallized_nuclei(nx, ny, i, j)) {
                    update_grid -> cell[nx][ny].grain_number = this -> cell[i][j].grain_number;
                    update_grid -> cell[nx][ny].dislocation_density = this -> cell[i][j].dislocation_density;
                    update_grid -> cell[nx][ny].N_recrystallized = 1;
                    update_grid -> cell[nx][ny].orientation = this -> cell[i][j].orientation;
                    propagated = 1;
                }
            }
        }
    }
    return propagated;
}

void deep_copy_grid(_grid_ *from, _grid_ *to) {
    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            to -> cell[i][j].orientation = from -> cell[i][j].orientation;
            to -> cell[i][j].N_recrystallized = from -> cell[i][j].N_recrystallized;
            to -> cell[i][j].grain_number = from -> cell[i][j].grain_number;
            to -> cell[i][j].gb_disp = from -> cell[i][j].gb_disp;
            to -> cell[i][j].dislocation_density = from ->cell[i][j].dislocation_density;
            to -> cell[i][j].gb_velocity = from -> cell[i][j].gb_velocity;
            to -> grain_num[i][j] = from -> grain_num[i][j];
        }
    }
    to -> grain_size_nc = from -> grain_size_nc;
    to -> v_max = from -> v_max;
    to -> p_max = from -> p_max;
    to -> p_avg = from -> p_avg;
}

void write_to_file(grain_cell array[GRID_SIZE][GRID_SIZE]) {
    // writing the file to be compatible with GNUPLOT (1st row and column are the axes values)
    // The indices have been inverted for the y-axis to match the origin and axes of the matrix with that of the plot.
    ofstream output("grid.dat");
    int i = 0, j = 0;
    ff(i, 0, GRID_SIZE + 1) {
        if (i == 0) {
            ff(k, 0, GRID_SIZE + 1) {
                output << k << " ";
            }
            output << "\n";
            continue;
        }
        ff(j, 0, GRID_SIZE + 1) {
            if (j == 0) {
                output << GRID_SIZE - i + 1 << " ";
                continue;
            }
            output << array[i - 1][j - 1].grain_number << " ";
        }
        output << "\n";
    }
    output.close();
}

void print_grid(grain_cell array[GRID_SIZE][GRID_SIZE]) {
    int i = 0, j = 0;
    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            cout << array[i][j].grain_number << " ";
        }
        cout << "\n";
    }
}

#endif