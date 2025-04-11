#include "main.h"

using namespace std;

int main()
{
    std::string path;

    SIMULATOR sim=SIMULATOR("mapas/ypacarai_map.csv");
    sim.reset();
    // sim.step();
    vector<string> msg {"Hello", "C++", "World"};

    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
}