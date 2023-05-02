#ifndef _DRX_GRID_
#define _DRX_GRID_

#include "helpers.h"
#include <vector>
#include <random>

using namespace std;

std::default_random_engine generator((unsigned)time(0));

std::normal_distribution<double> distribution (DIS_DEN_MEAN, DIS_DEN_STD_DEV);

class grain_cell {
    public:
    float orientation;
    int N_recrystallized;
    int grain_number; 
    float gb_disp;
    double dislocation_density;
    float gb_velocity;

    void update_dislocation_density(float delta_str);

    void nucleate();

    grain_cell() {
        this -> grain_number = rand() % STATES + 1;
        this -> orientation = ((float)rand() / (float)RAND_MAX) * (float)MAX_ORIENTATION;
        this -> gb_disp = 0;
        this -> N_recrystallized = 0;
        this -> gb_velocity = 0;
        this -> dislocation_density = distribution(generator);
    }
};

class _grid_ {
    public: 
    grain_cell cell[GRID_SIZE][GRID_SIZE];
    int grain_num[GRID_SIZE][GRID_SIZE];
    float v_max = 0;
    float p_max = 0;
    double p_avg = 0;
    vector<int> grain_size_nc;

    void monteCarloInit(grain_cell array[GRID_SIZE][GRID_SIZE]);
    void voronoiInit(grain_cell array[GRID_SIZE][GRID_SIZE]);
    float calculate_velocities();
    void average_p();
    void reset_grain_variables();
    void set_grain_numbers();
    void saturate_grain(int, int, int, int, bool);
    int check_neighbours(int, int);
    float calculate_cell_velocity(float, float, float);
    bool consume_recrystallized_nuclei(int, int, int, int);
    bool propagate_grain_boundary(int, int, int, _grid_*);
    bool potential_nucleus(int, int);

    _grid_() {
        encoder = 0;
        while ((int)GRID_SIZE / (int)pow(10, encoder) != 0) {
            encoder++;
        }
        this -> grain_size_nc = {0};
        monteCarloInit(this -> cell);
        set_grain_numbers();
        this -> v_max = calculate_velocities();
        average_p();
        std::cout << "Maximum Grain Velocity: " << this -> v_max << "\n";
    }

};

void deep_copy_grid(_grid_*, _grid_*);
void write_to_file(grain_cell array[GRID_SIZE][GRID_SIZE]);

#endif