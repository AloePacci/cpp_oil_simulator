#include "basic_simulator.h"
using namespace std;

// #define LOG_EVERYTHING_SIM

SIMULATOR::SIMULATOR(std::string _filepath, double _dt, double _kw, double _kc, double _gamma, double _flow, int _number_of_sources): dt(_dt),kw(_kw),kc(_kc),gamma(_gamma),flow(_flow),number_of_sources(_number_of_sources)
{
    mapa = new MAP(_filepath);
    
    source_points = Eigen::MatrixXi(2,number_of_sources);

    y = Eigen::VectorXi::LinSpaced(mapa->ncols, 0, mapa->ncols).rowwise().replicate(mapa->nrows);
    x = Eigen::RowVectorXi::LinSpaced(mapa->nrows, 0, mapa->nrows).colwise().replicate(mapa->ncols);
    // mapa->print();
    
}


void SIMULATOR::reset(){
    /* RESET THE ENV VARIABLES*/
    if(!init){
        init=true;
    }
    done=false;

    //Generate the source points
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, mapa->visitable.size()-1); // define the range
    for(int i=0;i<number_of_sources;i++){
        int index=distr(gen);
        source_points(0,i)=mapa->visitable[index].first;
        source_points(1,i)=mapa->visitable[index].second;
    }

    //generate wind speed
    wind_speed = Eigen::VectorXd::Random(2)*2-Eigen::VectorXd::Ones(2);
    contamination_position=source_points;


    //Random current vector field
    std::uniform_int_distribution<> distr_mapx(0,mapa->ncols-1); // define the range
    std::uniform_int_distribution<> distr_mapy(0,mapa->nrows-1); // define the range
    std::uniform_int_distribution<> distr2(3, 100); // define the range

    Eigen::MatrixXi x0 = Eigen::MatrixXi::Constant(mapa->ncols, mapa->nrows, distr_mapx(gen));
    Eigen::MatrixXi y0 = Eigen::MatrixXi::Constant(mapa->ncols, mapa->nrows, distr_mapy(gen));
    Eigen::MatrixXd aux1=(x-x0).cast<double>();
    Eigen::MatrixXd aux2=(y-y0).cast<double>();
    u=(aux1*(EIGEN_PI/distr2(gen))).array().sin()*((aux2*(M_PI/distr2(gen))).array().cos());
    v=(aux1*(EIGEN_PI/distr2(gen))).array().cos()*((aux2*(M_PI/distr2(gen))).array().sin());

    //density map
    density = Eigen::MatrixXi::Zero(mapa->nrows,mapa->ncols); //reset to zeros

    #ifdef LOG_EVERYTHING_SIM
        cout << "map shape" << mapa->ncols  << " " << mapa->nrows<< endl;
        cout << "source_points shape" <<source_points.rows()  << " " << source_points.cols()<< endl;
        cout << "contamination_position shape" <<contamination_position.rows()  << " " << contamination_position.cols()<< endl;
        cout << "u shape" <<u.rows()  << " " << u.cols()<< endl;
        cout << "v shape" <<v.rows()  << " " << v.cols()<< endl;
    #endif
}


void SIMULATOR::step(){
    assert(init); //"Environment not initiated!"
    //generate new particles
    for(int i=0;i<source_points.cols();i++){
        if(source_fuel>0){
            //compute components of the particle movement
            Eigen::VectorXd v_random = (Eigen::VectorXd::Random(2)*2-Eigen::VectorXd::Ones(2)) * gamma;
            Eigen::VectorXd v_wind = kw*wind_speed;
            Eigen::VectorXd aux(2);
            aux << v(source_points(1,i),source_points(0,i)), u(source_points(1,i),source_points(0,i));
            Eigen::VectorXd v_current = kc *aux; 
            //add new position to the list
            Eigen::VectorXi vnew= source_points.col(i) + (dt * (v_wind+ v_current) + v_random).cast<int>();
            std::pair<int, int> item(vnew(0), vnew(1));
            //if particle is not visitable, dont update pos
            if ( std::find(mapa->visitable.begin(), mapa->visitable.end(), item) == mapa->visitable.end() ){ 
                vnew= source_points.col(i);
            }
            contamination_position.conservativeResize(contamination_position.rows(), contamination_position.cols()+1);
            contamination_position.col(contamination_position.cols()-1)=vnew;
        }
    }

    //update particles positions
    for(int i=0;i<contamination_position.cols();i++){
        //compute components of the particle movement
        Eigen::VectorXd v_random = (Eigen::VectorXd::Random(2)*2-Eigen::VectorXd::Ones(2)) * gamma;
        Eigen::VectorXd v_wind = kw*wind_speed;
        Eigen::VectorXd aux(2);
        aux << v(contamination_position(1,i),contamination_position(0,i)), u(contamination_position(1,i),contamination_position(0,i));
        Eigen::VectorXd v_current = kc *aux; 

        //add new position to the list
        Eigen::VectorXi vnew= contamination_position.col(i) + (dt * (v_wind+ v_current) + v_random).cast<int>();
        std::pair<int, int> item(vnew(0), vnew(1));
        //if particle is visitable
        if ( std::find(mapa->visitable.begin(), mapa->visitable.end(), item) != mapa->visitable.end() ){ 
            //if density is below max
            if(density(vnew(0),vnew(1))<max_contamination_value)
            contamination_position.col(i)=vnew;
        }
        if(vnew(0)>mapa->ncols || vnew(0) <0|| vnew(1) < 0 || vnew(1)>mapa->nrows){ //if out of bounds, remove particle
            this->removeColumn(contamination_position,i);
        }
    }
    this->calculate_density();
}


void SIMULATOR::removeColumn(Eigen::MatrixXi &matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


void SIMULATOR::calculate_density(void){
    density = Eigen::MatrixXi::Zero(mapa->nrows,mapa->ncols); //reset to zeros

    for(int i=0; i<contamination_position.cols();i++){
        density(contamination_position(0,i), contamination_position(1,i))+=1;
    }
}