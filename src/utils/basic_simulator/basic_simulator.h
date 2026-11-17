#ifndef BASIC_SIMULATOR_CLASS_H
#define BASIC_SIMULATOR_CLASS_H

#include <iostream>
#include <vector>
#include <string>
#include "map.h"
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <algorithm>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX      // windows.h defines min/max macros that break std::min / std::max (MSVC error C2589)
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
class SIMULATOR
{
private:

    ////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////PARAMETERS/////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    double dt;       //time step
    double kw;      //weight of the windspeed
    double kc;        //weight of the current speed
    double gamma;     //weight of the random movement (brownian movement term)
    double flow;     //flow of oil - Number of particles generated per time step
    int number_of_sources=1;
    int max_contamination_value = 5;
    int source_fuel = 10000;
    bool apply_forces_at_origin=false;


    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////ENVIRONMENT PARAMETERS///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    bool done = false;
    int random_seed;
    int im0 = 0;
    int im1 = 0;
    float flow_remainder=0;


    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////ENVIRONMENT VARIABLES///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    float contamination_speed = 0;
    float spill_directions = 0;
    float particles_speeds = 0;
    bool init = false;   // was uninitialised: reading it in reset() was undefined behaviour (macOS CI: "Environment not initiated!")
    int step_counter =0;
    // initialising standard deviation to 1.0 
    double sigma = 1.0; 
    double GKernel[5][5];
    double max_kernel=0.0;


    MAP * mapa = NULL;
    void removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove);
    Eigen::VectorXd get_closest_neighbourgh(Eigen::VectorXd position);


    //for random
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen;

    ////////////////////////////////////////
    /////////Triangular movement
    ////////////////////////////////////////
    bool triangular;
    double triangular_magnitude = 1;
    double triangular_dilution = 1.2;

    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////// DYNAMIC EXTENSIONS (all OFF by default) //////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // v0.1.0: the three time-step corrections (diffusion ∝ sqrt(dt), emission per unit of time, no skipped particle) are
    // ON by default, with reference dt 0.3. Legacy (bit-identical to the PyPI 0.0.1 wheel for the same seed):
    // set_dynamics(diffusion_dt_ref=-1, emission_dt_ref=-1, fix_removal_skip=false). The optional dynamics (variable
    // wind, moving source, decay) are off by default and draw from their own generator (gen_dyn), never from `gen`.
    std::mt19937 gen_dyn;                 // seeded at reset from random_seed
    double wind_variability = 0.0;        // OU noise intensity sigma of the wind (sim units / sqrt(time)); 0 = constant wind
    double wind_relaxation = 1.0 / 30.0;  // OU mean-reversion rate theta towards the prevailing wind (1 / time)
    double wind_max = 1.4142135623730951; // cap on |wind| (the envelope of the original U(-1,1)^2 draw)
    double source_speed = 0.0;            // self-propelled source speed (cells / time); 0 = fixed source
    double source_turn = 0.0;             // heading diffusion of the moving source (rad / sqrt(time))
    double source_drift = 0.0;            // passive drift of the source with kw*wind + kc*current (0 = none, 1 = like oil)
    double decay_rate = 0.0;              // decay / evaporation hazard (1 / time): p_remove = 1 - exp(-rate * dt)
    double evaporable_fraction = 1.0;     // fraction of new particles that can decay (light fraction of the oil)
    double diffusion_dt_ref = 0.3;        // > 0: random walk scaled by sqrt(dt / diffusion_dt_ref) (Brownian); <= 0: legacy
    double emission_dt_ref = 0.3;         // > 0: emission per unit TIME (flow*dt/emission_dt_ref particles/step); <= 0: legacy per step
    bool fix_removal_skip = true;
    bool density_log = false;             // get_normalized_density: log1p(n)/log1p(max) instead of n/max (set_density)         // true: a particle removed at the border no longer makes the next one skip its step

    Eigen::MatrixXd source_pos;           // continuous source positions (2 x n), source_points = floor(source_pos)
    Eigen::VectorXd source_heading;       // heading of each moving source (rad)
    std::vector<char> evaporable;         // per particle, aligned with the columns of contamination_position

    Eigen::VectorXd draw2(std::uniform_real_distribution<> &d);  // 2 draws, first -> component 1 (PyPI wheel order)
    double diffusion_scale() const;
    void removeParticle(unsigned int i);  // removeColumn + keeps `evaporable` aligned
    void step_wind();
    void step_sources();
    void step_decay();
    bool visitable_cell(double r, double c) const;
public:
    Eigen::VectorXd wind_mean;            // prevailing wind (OU mean); = initial wind unless set_wind(.., true)
    int step_count() const { return step_counter; }
    int remaining_fuel() const { return source_fuel; }
    int particle_count() const { return (int)contamination_position.cols(); }

    void set_dynamics(double _wind_variability, double _wind_relaxation, double _wind_max, double _source_speed,
                      double _source_turn, double _source_drift, double _decay_rate, double _evaporable_fraction,
                      double _diffusion_dt_ref, double _emission_dt_ref, bool _fix_removal_skip);
    void warmup(int n_steps);             // run n_steps simulator steps (start the spill before the mission)
    void set_wind(Eigen::VectorXd w, bool also_mean = true);
    void set_source_fuel(int fuel);       // finite release (puff): emission stops when the fuel runs out
    void set_density(int max_per_cell, bool log_scale);   // crowding cap per cell + normalisation of the density map
    int get_max_per_cell() const { return max_contamination_value; }
    bool get_density_log() const { return density_log; }
    Eigen::MatrixXd get_source_position() const { return source_pos; }

    Eigen::MatrixXi source_points ;
    Eigen::MatrixXd contamination_position;
    Eigen::MatrixXi density;
    Eigen::MatrixXi x;
    Eigen::MatrixXi y;
    Eigen::MatrixXd v;
    Eigen::MatrixXd u;
    Eigen::VectorXd wind_speed;

    void reset(int _seed=-1, Eigen::MatrixXi _source_points_pos=Eigen::MatrixXi(), Eigen::MatrixXd _v=Eigen::MatrixXd(), Eigen::MatrixXd _u=Eigen::MatrixXd(), Eigen::VectorXd _wind_speed=Eigen::VectorXd(), Eigen::MatrixXd _contamination_position=Eigen::MatrixXd()); //these params are for initialization
    void step();
    Eigen::MatrixXd get_normalized_density(bool gaussian = true);

    SIMULATOR(std::string _filepath, double _dt=10, double _kw=0.5, double _kc=1, double _gamma=1, double _flow=1, int _number_of_sources=3, int _max_contamination_value=5, int _source_fuel = 1000, int _random_seed=-1, bool _triangular=false);
    SIMULATOR(Eigen::MatrixXi _base_matrix, double _dt=10, double _kw=0.5, double _kc=1, double _gamma=1, double _flow=1, int _number_of_sources=3, int _max_contamination_value=5, int _source_fuel = 1000, int _random_seed=-1, bool _triangular=false);
};


#endif