
//class
// #include "individual.hpp"
#include "ga.hpp"

//defines
#define CITY_SIZE 10
#define POP_SIZE 5

//variables
struct Point {
    float x, y;
};

// Los indices 0,1,3... son los genes
struct Point map[CITY_SIZE];