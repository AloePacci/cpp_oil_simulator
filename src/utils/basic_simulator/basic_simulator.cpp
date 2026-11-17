#include "basic_simulator.h"
#include <stdexcept>
using namespace std;

// #define LOG_EVERYTHING_SIM

SIMULATOR::SIMULATOR(std::string _filepath, double _dt, double _kw, double _kc, double _gamma, double _flow, int _number_of_sources, int _max_contamination_value, int _source_fuel, int _random_seed, bool _triangular): dt(_dt),kw(_kw),kc(_kc),gamma(_gamma),flow(_flow),number_of_sources(_number_of_sources), max_contamination_value(_max_contamination_value), source_fuel(_source_fuel), random_seed(_random_seed), triangular(_triangular)
{
    mapa = new MAP(_filepath);
    source_points = Eigen::MatrixXi(2,number_of_sources);
    y = Eigen::VectorXi::LinSpaced(mapa->ncols, 0, mapa->ncols).rowwise().replicate(mapa->nrows);
    x = Eigen::RowVectorXi::LinSpaced(mapa->nrows, 0, mapa->nrows).colwise().replicate(mapa->ncols);
    // mapa->print();

    if(random_seed == -1){
        gen =  std::mt19937(rd());
        cout << "no seed specified (-1), using random seed" << endl;
    }else{
        gen =  std::mt19937(random_seed);
    }

    #ifdef LOG_EVERYTHING_SIM
        cout << "found map of " << mapa->ncols << " colums and " << mapa->nrows << " rows" << endl;
        mapa->print();
    #endif


    // generating 5x5 kernel 
    double r, s = 2.0 * sigma * sigma; 
    for (int x = -2; x <= 2; x++) { 
        for (int y = -2; y <= 2; y++) { 
            r = sqrt(x * x + y * y); 
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s); 
            max_kernel += GKernel[x + 2][y + 2]; 
        } 
    } 
}

SIMULATOR::SIMULATOR(Eigen::MatrixXi _base_matrix, double _dt, double _kw, double _kc, double _gamma, double _flow, int _number_of_sources, int _max_contamination_value, int _source_fuel, int _random_seed, bool _triangular): dt(_dt),kw(_kw),kc(_kc),gamma(_gamma),flow(_flow),number_of_sources(_number_of_sources), max_contamination_value(_max_contamination_value), source_fuel(_source_fuel), random_seed(_random_seed), triangular(_triangular)
{
    mapa = new MAP(_base_matrix);
    source_points = Eigen::MatrixXi(2,number_of_sources);
    y = Eigen::VectorXi::LinSpaced(mapa->ncols, 0, mapa->ncols).rowwise().replicate(mapa->nrows);
    x = Eigen::RowVectorXi::LinSpaced(mapa->nrows, 0, mapa->nrows).colwise().replicate(mapa->ncols);
    // mapa->print();

    if(random_seed == -1){
        gen =  std::mt19937(rd());
        cout << "no seed specified (-1), using random seed" << endl;
    }else{
        gen =  std::mt19937(random_seed);
    }

    #ifdef LOG_EVERYTHING_SIM
        cout << "found map of " << mapa->ncols << " colums and " << mapa->nrows << " rows" << endl;
        mapa->print();
    #endif


    // generating 5x5 kernel 
    double r, s = 2.0 * sigma * sigma; 
    for (int x = -2; x <= 2; x++) { 
        for (int y = -2; y <= 2; y++) { 
            r = sqrt(x * x + y * y); 
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s); 
            max_kernel += GKernel[x + 2][y + 2]; 
        } 
    } 
}


void SIMULATOR::reset(int _seed, Eigen::MatrixXi _source_points_pos, Eigen::MatrixXd _v, Eigen::MatrixXd _u, Eigen::VectorXd _wind_speed, Eigen::MatrixXd _contamination_position){ 
    if(_seed != -1){
        random_seed = _seed;
        gen =  std::mt19937(random_seed);
    #ifdef LOG_EVERYTHING_SIM
        cout << "Reset with seed " <<random_seed  << endl;
    }else{
        cout << "Reset with same seed " <<random_seed  << endl;
    #endif
    }
    /* RESET THE ENV VARIABLES*/
    init = true;

    done=false;
    step_counter =0;

    //Generate the source points
    std::uniform_int_distribution<> distr(0, mapa->visitable.size()-1); // define the range
    
    if (_source_points_pos.size()==0){
        for(int i=0;i<number_of_sources;i++){
            int index=distr(gen);
            source_points(0,i)=mapa->visitable[index].first;
            source_points(1,i)=mapa->visitable[index].second;
        }
    }else{
        if(_source_points_pos.cols()!=number_of_sources || _source_points_pos.rows()!=2){
            char message[100];
            sprintf(message, "bad source points shape %dx%d when it should be %dx2", _source_points_pos.cols(), _source_points_pos.rows(), number_of_sources);
            cout << message << endl;
            throw std::runtime_error(message);
        }
        for (int i=0;i<_source_points_pos.cols();i++){
            bool valid=false;
            for(auto point: mapa->visitable){
                if(point.first > mapa->nrows-1 || point.second > mapa->ncols-1){
                    char message[100];
                    sprintf(message, "bad source point position %d,%d is out of bounds", point.first, point.second);
                    cout << message << endl;
                    throw std::runtime_error(message);
                }
                if(point.first == _source_points_pos(0,i) && point.second == _source_points_pos(1,i)){
                    valid=true;
                    break;
                }
            }
            if(!valid){
                char message[100];
                sprintf(message, "bad source point position %d,%d is not visitable", _source_points_pos(0,i), _source_points_pos(1,i));
                cout << message << endl;
                throw std::runtime_error(message);
            }
        }
        source_points = _source_points_pos;
    }

    #ifdef LOG_EVERYTHING_SIM
        cout << "source_points shape" <<source_points.rows()  << "x" << source_points.cols()<< endl;
        cout << "source_points: " <<source_points<< endl;
    #endif

    std::uniform_real_distribution<> dist_wind(-1, 1); // define the range
    
    //generate wind speed
    if(_wind_speed.size()==0){
        wind_speed = Eigen::VectorXd(2);
        {   // explicit draw order (first draw -> component 1), as the PyPI 0.0.1 wheel: see draw2()
            double w1 = dist_wind(gen), w0 = dist_wind(gen);
            wind_speed << w0, w1;
        }
    }else if(_wind_speed.size()==2){
        wind_speed = _wind_speed;
    }else{
        char message[100];
        sprintf(message, "bad wind speed shape %dx%d when it should be 2x1", _wind_speed.cols(), _wind_speed.rows());
        cout << message << endl;
        throw std::runtime_error(message);
    }

    #ifdef LOG_EVERYTHING_SIM
        cout << "wind speed is [" << dt * kw *(wind_speed.transpose()) << "]" <<endl;
    #endif

    //initialize
    if (_contamination_position.size()!=0){
        if(_contamination_position.rows()!=2){
            char message[100];
            sprintf(message, "bad contamination position shape %dx%d when it should be 2xN", _contamination_position.cols(), _contamination_position.rows());
            cout << message << endl;
            throw std::runtime_error(message);
        }
        contamination_position=_contamination_position;
    }else{
        contamination_position=source_points.cast<double>();
    }

    // dynamic extensions: own generator (does not touch `gen`), prevailing wind, continuous source state
    gen_dyn = (random_seed == -1) ? std::mt19937(rd()) : std::mt19937((unsigned int)random_seed ^ 0x9E3779B9u);
    wind_mean = wind_speed;
    source_pos = source_points.cast<double>();
    source_heading = Eigen::VectorXd(number_of_sources);
    {
        std::uniform_real_distribution<> dist_head(-M_PI, M_PI);
        for (int i = 0; i < number_of_sources; i++) source_heading(i) = dist_head(gen_dyn);
    }
    evaporable.assign(contamination_position.cols(), 1);    // also for user-given contamination positions

    if( _v.size()!=0 || _u.size()!=0){
        if(_u.cols()!=mapa->ncols || _u.rows()!=mapa->nrows || _v.cols()!=mapa->ncols || _v.rows()!=mapa->nrows)
        {
            char message[100];
            sprintf(message, "bad u or v shape %dx%d  %dx%d when it should be %dx%d", _u.cols(), _u.rows(), _v.cols(), _v.rows(), mapa->ncols, mapa->nrows);
            cout << message << endl;
            throw std::runtime_error(message);
        }
        u=_u;
        v=_v;
    }else if(!triangular){
        //Random current vector field
        std::uniform_int_distribution<> distr_mapx(0,mapa->ncols-1); // define the range
        std::uniform_int_distribution<> distr_mapy(0,mapa->nrows-1); // define the range
        std::uniform_real_distribution<> distr2(3, 100); // define the range

        Eigen::MatrixXi x0 = Eigen::MatrixXi::Constant(mapa->ncols, mapa->nrows, distr_mapx(gen));
        Eigen::MatrixXi y0 = Eigen::MatrixXi::Constant(mapa->ncols, mapa->nrows, distr_mapy(gen));
        Eigen::MatrixXd aux1=(x-x0).cast<double>();
        Eigen::MatrixXd aux2=(y-y0).cast<double>();
        u=(aux1*(EIGEN_PI/distr2(gen))).array().sin()*((aux2*(M_PI/distr2(gen))).array().cos());
        v=-(aux1*(EIGEN_PI/distr2(gen))).array().cos()*((aux2*(M_PI/distr2(gen))).array().sin());
    }else{ //currents generated from start point
        u = Eigen::MatrixXd::Constant(mapa->ncols, mapa->nrows,0);
        v = Eigen::MatrixXd::Constant(mapa->ncols, mapa->nrows,0);
        double xt= source_points(0,0);
        double yt= source_points(1,0);
        for(int m=1;m<mapa->ncols-1;m++){
            for(int n=1;n<mapa->nrows-1;n++){
                if(n-xt == 0 && m-yt ==0){
                    std::uniform_int_distribution<> distr_start(-1, 1); // define the range
                    u(m,n)=distr_start(gen)*triangular_magnitude;
                    v(m,n)=distr_start(gen)*triangular_magnitude;
                }else{
                    if(abs(xt-n)>abs(m-yt)){
                        v(m,n)=triangular_magnitude*(n-xt)/pow(abs(n-xt),triangular_dilution);
                    }else{
                        u(m,n)=triangular_magnitude*(m-yt)/pow(abs(m-yt),triangular_dilution);
                    }
                }
            }
        }
    }

    //density map
    density = Eigen::MatrixXi::Constant(mapa->nrows,mapa->ncols,0); //initialize to -1
    //put to 0 map points
    // for (auto point: mapa->visitable){
    //     density(point.first,point.second)=0;
    // }
    // initialize first particles in density
    for(int i=0; i<contamination_position.cols();i++){
        density((int)contamination_position(0,i),(int)contamination_position(1,i))++;
    }


    #ifdef LOG_EVERYTHING_SIM
        cout << "map shape" << mapa->ncols  << " " << mapa->nrows<< endl;
        cout << "contamination_position shape" <<contamination_position.rows()  << " " << contamination_position.cols()<< endl;
        cout << "u shape" <<u.rows()  << " " << u.cols()<< endl;
        cout << "v shape" <<v.rows()  << " " << v.cols()<< endl;
        // cout << "Visitable:" << endl;
        // for(auto i: mapa->visitable){
        //     cout << "[" << i.first << "," << i.second << "], ";
        // } 
        cout << endl;
    #endif
}


void SIMULATOR::step(){
    if(!init) {
        throw std::runtime_error("Environment not initiated!");
    }
    //update particles positions
    std::uniform_real_distribution<> distrandom(-2, 2); // define the range
    for(int i=0;i<contamination_position.cols();i++){
        //compute components of the particle movement
        Eigen::VectorXd v_random = Eigen::VectorXd(2);
        v_random = draw2(distrandom) * (gamma * diffusion_scale());
        // else{
        //     v_random = _random.row(i)*gamma;
        // }
        Eigen::VectorXd v_wind = kw*wind_speed;
        Eigen::VectorXd aux(2);
        aux << v((int)contamination_position(1,i),(int)contamination_position(0,i)), u((int)contamination_position(1,i),(int)contamination_position(0,i));
        Eigen::VectorXd v_current = kc *aux; 

        //add new position to the list
        Eigen::VectorXd vnew= contamination_position.col(i) + (dt * (v_wind+ v_current) + v_random);
        #ifdef LOG_EVERYTHING_SIM
        cout << "wind [" << v_wind.transpose() << "]    update position[" << (dt * (v_wind+ v_current) + v_random).cast<int>().transpose() << "]    current [" << v_current.transpose() << "]     random["  <<  v_random.transpose() <<"]" <<endl;
        #endif
        //if particle is not visitable
        if((int)vnew(0)>mapa->nrows-1 || (int)vnew(0) <0|| (int)vnew(1) < 0 || (int)vnew(1)>mapa->ncols-1){ //if out of bounds, remove particle
            #ifdef LOG_EVERYTHING_SIM
            cout << "we got to remove " << contamination_position.col(i).transpose() << endl;
            #endif
            density((int)contamination_position(0,i),(int)contamination_position(1,i))--;
            if(density((int)contamination_position(0,i),(int)contamination_position(1,i)) <0){
                cout << "something weird happened, erased particle where should be none, move out at " << step_counter  << contamination_position.col(i).transpose() << endl;
            }
            this->removeParticle(i);
            if (fix_removal_skip) { i--; continue; }   // the next particle moved into column i: process it too
        }
        else if ( mapa->map.at((int)vnew(1)).second[(int)vnew(0)] ==1){ //if position is not visitable, position wont be updated
            density((int)contamination_position(0,i),(int)contamination_position(1,i))--; //update density and position
            if(density((int)contamination_position(0,i),(int)contamination_position(1,i)) <0){
                cout << "something weird happened, erased particle where should be none, move in at " << step_counter  <<"  " << contamination_position.col(i).transpose() << endl;
                density((int)contamination_position(0,i),(int)contamination_position(1,i))=0;
            }
            
            if(density((int)(vnew(0)),(int)vnew(1))>=max_contamination_value){ //check density before creating new particle
                vnew= get_closest_neighbourgh(vnew);//look for inmediate position next to it if it is full
            }
            contamination_position.col(i)=vnew;
            density((int)contamination_position(0,i),(int)contamination_position(1,i))++;
        }
    }
    if (decay_rate > 0) step_decay();        // decay / evaporation of the existing particles
    //generate new particles
    for(int i=0;i<source_points.cols();i++){
        if(source_fuel>0){
            // original: ~flow particles per STEP (whatever dt). emission_dt_ref > 0: flow is a rate per unit of TIME,
            // particles per step = flow * dt / emission_dt_ref (same mean rate as the original at dt = emission_dt_ref)
            flow_remainder += (emission_dt_ref > 0) ? (float)(flow * dt / emission_dt_ref) : (float)flow;
            while(((emission_dt_ref > 0) ? (flow_remainder >= 1) : (flow_remainder*dt>1)) && source_fuel>0){ //generate as many particles as flow says
                Eigen::VectorXd v_random = Eigen::VectorXd(2);
                v_random = draw2(distrandom) * (gamma * diffusion_scale());
                // else{
                //     v_random = _random.row(i+contamination_position.cols())*gamma;
                // }
                //compute components of the particle movement
                Eigen::VectorXd v_wind = kw*wind_speed;
                Eigen::VectorXd aux(2);
                aux << v(source_points(1,i),source_points(0,i)), u(source_points(1,i),source_points(0,i));
                Eigen::VectorXd v_current = kc *aux; 
                //add new position to the list
                Eigen::VectorXd vnew(2);
                if (apply_forces_at_origin){
                    vnew= source_points.col(i).cast<double>() + (dt * (v_wind+ v_current) + v_random);
                }else{
                    vnew= source_points.col(i).cast<double>() +  v_random;
                }
                // std::pair<int, int> item(vnew(0), vnew(1));
                //if particle is not visitable, dont update pos
                if((int)vnew(0)>mapa->nrows-1 || (int)vnew(0) <0|| (int)vnew(1) < 0 || (int)vnew(1)>mapa->ncols-1){
                    cout << "something really bad happened with simulator, particle spawned outside map" << endl;
                    flow_remainder--; //to avoid infinite wait. if there is no fuel to generate, dont
                    continue;
                }
                if ( mapa->map.at((int)vnew(1)).second[(int)vnew(0)] ==0){ 
                    vnew= source_points.col(i).cast<double>();
                }
                contamination_position.conservativeResize(contamination_position.rows(), contamination_position.cols()+1);
                if(density((int)(vnew(0)),(int)vnew(1))>=max_contamination_value){ //check density before creating new particle
                    vnew= get_closest_neighbourgh(vnew);//look for inmediate position next to it if it is full
                }
                contamination_position.col(contamination_position.cols()-1)=vnew;
                {   // evaporable flag of the new particle (extra draw only when decay with a residual fraction is on)
                    std::uniform_real_distribution<> u01(0.0, 1.0);
                    bool flag = (decay_rate > 0 && evaporable_fraction < 1.0) ? (u01(gen_dyn) < evaporable_fraction) : true;
                    evaporable.push_back(flag ? 1 : 0);
                }
                density((int)contamination_position(0,contamination_position.cols()-1),(int)contamination_position(1,contamination_position.cols()-1))++;
                source_fuel--;
                flow_remainder--;
            }
        }
    }
    if (source_speed > 0 || source_drift > 0) step_sources();   // move the plume origin for the next step
    if (wind_variability > 0) step_wind();                      // wind for the next step
    step_counter++;
}


void SIMULATOR::removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


Eigen::MatrixXd SIMULATOR::get_normalized_density(bool gaussian){
    Eigen::MatrixXd density_d = Eigen::MatrixXd::Zero(mapa->nrows,mapa->ncols);
    if(gaussian==false){
        for (auto point: mapa->visitable){
            density_d(point.first,point.second)=((double)density(point.first,point.second))/((double)(max_contamination_value));
        }
    }else{
        for (auto point: mapa->visitable){
            double kernel=max_kernel;
            for (int x = -2; x <= 2; x++) { 
                for (int y = -2; y <= 2; y++) { 
                    if(point.first+x >1 && point.first+x < mapa->nrows-2 && point.second+y >1 && point.second+y <mapa->ncols -1){
                        if (mapa->map.at(point.second+y).second[point.first+x] ==0)
                            kernel-=GKernel[x+2][y+2];
                        else
                            density_d(point.first,point.second)+=((double)density(point.first+x,point.second+y))*GKernel[x+2][y+2]/((double)(max_contamination_value));
                    }else{
                        kernel-=GKernel[x+2][y+2];
                    }
                }
            }
            density_d(point.first,point.second)/=kernel;
        }
    }
    if (density_log){   // logarithmic scale: log(1 + n) / log(1 + max), n = (smoothed) particles per cell
        double lm = std::log1p((double)max_contamination_value);
        for (auto point: mapa->visitable){
            double n = density_d(point.first,point.second) * (double)max_contamination_value;
            density_d(point.first,point.second) = std::log1p(n > 0.0 ? n : 0.0) / lm;
        }
    }
    return density_d;
}


Eigen::VectorXd SIMULATOR::get_closest_neighbourgh(Eigen::VectorXd position){
    Eigen::VectorXd test(2);
    int radius=1;
    while(true){ //get contouring points
        if(radius>mapa->ncols/3+mapa->nrows/3){
            cout << "looked for too long" <<endl;
            return position;
        }
        Eigen::MatrixXd available = Eigen::MatrixXd::Constant(2, 1, -1);;
        for(int i=-radius;i<=radius;i++){
            for(int j=-radius;j<=radius;j++){
                if(int(sqrt(i*i+j*j))==radius){
                    if((int)(position(0)+i)>mapa->nrows-1 || (int)(position(0)+i) <0|| (int)(position(1)+j) < 0 || (int)(position(1)+j)>mapa->ncols-1)
                        continue; //if not visitable, look for next
                    if(mapa->map.at((int)(position(1)+j)).second[(int)(position(0)+i)] !=1) //if it is not visitable
                        continue;
                    if(density((int)(position(0)+i),(int)(position(1)+j))<max_contamination_value){//is  a valid position
                        available(0,available.cols()-1)=position(0)+i;
                        available(1,available.cols()-1)=position(1)+j;
                        available.conservativeResize(2, available.cols()+1);
                    }
                }
            }
        }
        radius++;
        if(available.cols()==1 && available(0,0) == -1){ //if we didnt find any available point
            continue;
        }
        //we found a point, pick a random position
        std::uniform_int_distribution<> distr_available(0,available.cols()-2);
        return available.col(distr_available(gen));
    }
    return position;
}

//TODO: check integrity of data

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// DYNAMIC EXTENSIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd SIMULATOR::draw2(std::uniform_real_distribution<> &d){
    // The PyPI 0.0.1 wheel evaluates `v << d(gen), d(gen)` right to left (GCC 10, pre-C++17 rules): the FIRST draw
    // lands in component 1. Written explicitly so that the result does not depend on the compiler / standard.
    double r1 = d(gen);
    double r0 = d(gen);
    Eigen::VectorXd out(2);
    out << r0, r1;
    return out;
}

double SIMULATOR::diffusion_scale() const {
    // Brownian motion: displacement std grows with sqrt(dt). dt_ref = the dt at which gamma was tuned.
    return diffusion_dt_ref > 0 ? std::sqrt(dt / diffusion_dt_ref) : 1.0;
}

void SIMULATOR::removeParticle(unsigned int i){
    removeColumn(contamination_position, i);
    if (i < evaporable.size()) evaporable.erase(evaporable.begin() + i);
}

bool SIMULATOR::visitable_cell(double r, double c) const {
    if ((int)r < 0 || (int)c < 0 || (int)r > mapa->nrows - 1 || (int)c > mapa->ncols - 1 || r < 0 || c < 0) return false;
    return mapa->map.at((int)c).second[(int)r] == 1;
}

void SIMULATOR::set_dynamics(double _wind_variability, double _wind_relaxation, double _wind_max, double _source_speed,
                             double _source_turn, double _source_drift, double _decay_rate, double _evaporable_fraction,
                             double _diffusion_dt_ref, double _emission_dt_ref, bool _fix_removal_skip){
    if (_wind_variability < 0 || _wind_relaxation < 0 || _wind_max <= 0 || _source_speed < 0 || _source_turn < 0 ||
        _source_drift < 0 || _decay_rate < 0 || _evaporable_fraction < 0 || _evaporable_fraction > 1)
        throw std::runtime_error("set_dynamics: negative rate/speed, wind_max <= 0 or evaporable_fraction outside [0, 1]");
    wind_variability = _wind_variability;
    wind_relaxation = _wind_relaxation;
    wind_max = _wind_max;
    source_speed = _source_speed;
    source_turn = _source_turn;
    source_drift = _source_drift;
    decay_rate = _decay_rate;
    evaporable_fraction = _evaporable_fraction;
    diffusion_dt_ref = _diffusion_dt_ref;
    emission_dt_ref = _emission_dt_ref;
    fix_removal_skip = _fix_removal_skip;
}

void SIMULATOR::warmup(int n_steps){
    if (!init) throw std::runtime_error("Environment not initiated! call reset() before warmup()");
    for (int k = 0; k < n_steps; k++) step();
}

void SIMULATOR::set_wind(Eigen::VectorXd w, bool also_mean){
    if (w.size() != 2) throw std::runtime_error("set_wind: wind must have 2 components");
    wind_speed = w;
    if (also_mean) wind_mean = w;
}

void SIMULATOR::set_source_fuel(int fuel){
    if (fuel < 0) throw std::runtime_error("set_source_fuel: fuel must be >= 0");
    source_fuel = fuel;
}

void SIMULATOR::set_density(int max_per_cell, bool log_scale){
    // max_per_cell: crowding cap (a full cell pushes newcomers to a free neighbour) AND the value that maps to 1.
    // A larger cap keeps the thickness of the slick (less saturation); log_scale spreads the low concentrations.
    if (max_per_cell < 1) throw std::runtime_error("set_density: max_per_cell must be >= 1");
    max_contamination_value = max_per_cell;
    density_log = log_scale;
}

void SIMULATOR::step_wind(){
    // Ornstein-Uhlenbeck around the prevailing wind (Euler-Maruyama), magnitude capped at wind_max.
    // Stationary std per component = wind_variability / sqrt(2 * wind_relaxation); correlation time = 1 / wind_relaxation.
    std::normal_distribution<> n01(0.0, 1.0);
    double sq = std::sqrt(dt);
    for (int k = 0; k < 2; k++)
        wind_speed(k) += wind_relaxation * (wind_mean(k) - wind_speed(k)) * dt + wind_variability * sq * n01(gen_dyn);
    double m = wind_speed.norm();
    if (m > wind_max) wind_speed *= wind_max / m;
}

void SIMULATOR::step_sources(){
    // Moving plume origin: self-propelled (speed + heading random walk, bounces on land / border)
    // plus passive drift with the same forcing as the oil (kw * wind + kc * current).
    std::normal_distribution<> n01(0.0, 1.0);
    std::uniform_real_distribution<> bounce(-M_PI / 4, M_PI / 4);
    for (int i = 0; i < source_pos.cols(); i++){
        Eigen::Vector2d d(0.0, 0.0);
        if (source_speed > 0){
            source_heading(i) += source_turn * std::sqrt(dt) * n01(gen_dyn);
            d += source_speed * dt * Eigen::Vector2d(std::cos(source_heading(i)), std::sin(source_heading(i)));
        }
        if (source_drift > 0){
            int r = (int)source_pos(0, i), c = (int)source_pos(1, i);
            Eigen::Vector2d forcing(kw * wind_speed(0) + kc * v(c, r), kw * wind_speed(1) + kc * u(c, r));
            d += source_drift * dt * forcing;
        }
        Eigen::Vector2d p = source_pos.col(i) + d;
        if (visitable_cell(p(0), p(1))){
            source_pos.col(i) = p;
        } else {
            source_heading(i) += M_PI + bounce(gen_dyn);   // turn back (stays in place this step)
        }
        source_points(0, i) = (int)source_pos(0, i);
        source_points(1, i) = (int)source_pos(1, i);
    }
}

void SIMULATOR::step_decay(){
    // Each evaporable particle disappears with probability 1 - exp(-decay_rate * dt) per step (exponential mass loss,
    // half-life ln2 / decay_rate time units; the non-evaporable fraction is a persistent residue).
    std::uniform_real_distribution<> u01(0.0, 1.0);
    double p = 1.0 - std::exp(-decay_rate * dt);
    for (int i = (int)contamination_position.cols() - 1; i >= 0; i--){
        if ((size_t)i < evaporable.size() && !evaporable[i]) continue;
        if (u01(gen_dyn) < p){
            int r = (int)contamination_position(0, i), c = (int)contamination_position(1, i);
            if (density(r, c) > 0) density(r, c)--;
            removeParticle(i);
        }
    }
}
