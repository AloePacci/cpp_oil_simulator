#ifndef MAP_H
#define MAP_H
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <Eigen/Dense>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
class MAP
{
private:
    std::string mappath;
    std::string mapname;
    bool exists();
    void read_csv();
    char delimiter=',';

public:
    MAP(std::string _filepath);
    MAP(Eigen::MatrixXi _map);
    std::vector<std::pair<int, std::vector<int>>> map;
    void print() const;
    int ncols;
    int nrows;

    std::vector<std::pair<int, int>> visitable;
};

#endif