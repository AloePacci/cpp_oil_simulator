#ifndef BASIC_SIMULATOR_CLASS_H
#define BASIC_SIMULATOR_CLASS_H

#include <iostream>
#include <vector>
#include <string>
#include "map.h"
#include <eigen3/Eigen/Dense>
#include <random>
#include <unistd.h>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

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
    bool done;
    int random_seed;
    int im0;
    int im1;
    float flow_remainder=0;


    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////ENVIRONMENT VARIABLES///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    float contamination_speed;
    float spill_directions;
    float particles_speeds;
    bool init;
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
public:
    Eigen::MatrixXi source_points ;
    Eigen::MatrixXd contamination_position;
    Eigen::MatrixXi density;
    Eigen::MatrixXi x;
    Eigen::MatrixXi y;
    Eigen::MatrixXd v;
    Eigen::MatrixXd u;
    Eigen::VectorXd wind_speed;

    void reset(int _seed=-1); //these params are for initialization
    void step();
    Eigen::MatrixXd get_normalized_density(bool gaussian = true);

    SIMULATOR(std::string _filepath, double _dt=10, double _kw=0.5, double _kc=1, double _gamma=1, double _flow=1, int _number_of_sources=3, int _max_contamination_value=5, int _source_fuel = 1000, int _random_seed=-1, bool _triangular=false);
    SIMULATOR(Eigen::MatrixXi _base_matrix, double _dt=10, double _kw=0.5, double _kc=1, double _gamma=1, double _flow=1, int _number_of_sources=3, int _max_contamination_value=5, int _source_fuel = 1000, int _random_seed=-1, bool _triangular=false);
};


#endif