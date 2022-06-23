#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cmath>
using std::sqrt;
using std::pow;

#include <fstream>
using std::ifstream;

#include <numeric>
#include <algorithm>

#include "main.hpp"

// struct Point map[CITY_SIZE];
vector<Point> map;

/**************** UTILS ********************/

void readFile() {
    string line;
    ifstream myfile ("data/dj38.tsp");

    int value;
    if (myfile.is_open()) {
        while ( getline (myfile,line) ) {
            // cout << line.substr(0,9) <<endl;
            if(line.substr(0,9) == "DIMENSION") {
                value = std::stoi(line.substr(12));
                cout << "value:" << value << endl;
                continue;
            }
            if (line.substr(0,18) == "NODE_COORD_SECTION")
                break;
        }

        while ( getline (myfile,line) && line[0] != 'E') {
            string s = line; 
            string delimiter = " ";
            size_t pos_start = 0, pos_end, delim_len = delimiter.length();
            string token;
            vector<string> res;

            while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
                token = s.substr (pos_start, pos_end - pos_start);
                pos_start = pos_end + delim_len;
                res.push_back (token);
            }

            res.push_back (s.substr (pos_start));
            // for (int i = 0; i < value; i++) {
            //     Point p;
            //     p.x = std::stof(res[1]);
            //     p.y = std::stof(res[2]);
            //     map.push_back(p);
            // }
            Point p;
            p.x = std::stof(res[1]);
            p.y = std::stof(res[2]);
            map.push_back(p);
            // cout << line << endl;
        }

        // for (size_t i = 0; i < value; i++) {
        //     cout << i <<") x: " << map[i].x << " y: " << map[i].y << endl;
        // }

        myfile.close();
    }
    else cout << "Unable to open file";
}

size_t randNumber (size_t start, size_t end) {
	size_t range = (end-start)+1;
	size_t random_int = start+(rand()%range);
	return random_int;
}


void burbleSort(vector<Individual>& population) {
    for (size_t i = 0; i < population.size(); i++) {
        for (size_t j = 0; j < population.size(); j++) {
            if (population[i].getFitness() < population[j].getFitness()) {
                Individual temp = population[i];
                population[i] = population[j];
                population[j] = temp;
            }
        }        
    }    
}

bool isValueinRange(vector<size_t> vec, size_t value, size_t start, size_t end) {
    for (size_t i = start; i <= end; i++)
        if (vec[i] == value)
            return true;

    return false;
}

int getBestFitnessIndex(const  vector<Individual>& population) {
    int best_index = 0;
    for (size_t i = 1; i < population.size()-1; i++) {
        if (population[best_index].getFitness() > population[i].getFitness())
            best_index = i;
    }
    return best_index;
}

Individual tournamentSelection(const vector<Individual>& population, int k, int range) {
    if (k > POP_SIZE) 
        k = POP_SIZE;

    int select[k];
    for (int i = 0; i < k; i++)
        select[i] = randNumber(0, range - 1);

    int bestFitness = select[0];
    for (int i = 1; i < k-1; i++) {
        //el menor fitness es el mejor
        if(population[bestFitness].getFitness() > population[select[i]].getFitness())
            bestFitness = select[i];
    }
    
    return population[bestFitness];
}

void singlePointCrossovers(const Individual& parent1, const Individual& parent2) {

}

void orderCrossover_OX(Individual& parent1, Individual& parent2) {
    size_t i_start = 0;
    size_t i_end = 0;
    
    do {
        i_start = randNumber(0, map.size() - 1);
        i_end   = randNumber(0, map.size() - 1);
    }while(i_start == i_end);
    
    if (i_start > i_end) {
        const int temp = i_start;
        i_start = i_end;
        i_end = temp;
    }

    // los nuevos genes con el cruzamiento
    vector<size_t> new_genes1;
    vector<size_t> new_genes2;
    new_genes1.resize(map.size());
    new_genes2.resize(map.size());

    // copiamos la seccion del padre 1
    for (size_t i = i_start; i <= i_end; i++)
        new_genes1[i] = parent1._genes[i];
    // copiamos la seccion del padre 2
    for (size_t i = i_start; i <= i_end; i++) 
        new_genes2[i] = parent2._genes[i];

    // se crea el hijo del padre1 con los genes del padre 2
    int index = 0;
    for (size_t i = 0; i < map.size(); i++)  {
        if (i >= i_start && i <= i_end) {
            new_genes1[i] = parent1._genes[i];
            continue;
        }
        if (!isValueinRange(new_genes1, parent2._genes[index], i_start, i_end)) {
            new_genes1[i] = parent2._genes[index];
            index++;
        } 
        else {
            index++;
            i--;
        }

    }

    // se crea el hijo del padre2 con los genes del padre 1
    index = 0;
    for (size_t i = 0; i < map.size(); i++)  {
        if (i >= i_start && i <= i_end) {
            new_genes2[i] = parent2._genes[i];
            continue;
        }
        if (!isValueinRange(new_genes2, parent1._genes[index], i_start, i_end)) {
            new_genes2[i] = parent1._genes[index];
            index++;
        } 
        else {
            index++;
            i--;
        }

    }
    // los nuevos hijos seran los nuevos genes            
    parent1._genes = new_genes1;
    parent2._genes = new_genes2;
}

// Individual pmxCrossovers(const Individual& parent1, const Individual& parent2) {
//     int p1 = randNumber(0, CITY_SIZE - 1);
//     int p2 = randNumber(0, CITY_SIZE - 1);
//     int pm = (p1 + p2) / 2;
// }

/**************** CLASS INDIVIDUAL ********************/
Individual::Individual() {
    _genes.resize(map.size());
    std::iota(_genes.begin(), _genes.end(), 0);
    initRandGenes();
}

Individual::Individual(const vector<size_t>& genes) {
    _genes.resize(map.size());
    for (size_t i = 0; i < _genes.size(); i++)
        _genes[i] = genes[i];
}

Individual::~Individual() { }

void Individual::initRandGenes() {
    for (size_t i = 0; i < map.size(); i++)
    {
        size_t rand_gene = randNumber(0, map.size() - 1);        
        size_t gene_temp = _genes[i];
        _genes[i] = _genes[rand_gene];
        _genes[rand_gene] = gene_temp;
    }
    // _fitness = calcRouteDistance();       
}

void Individual::calcFitness() {
    _fitness = calcRouteDistance();
}

void Individual::swapMutation() {
   // Mutamos el individuo intercambiando dos de sus genes de
   // forma aleatoria

    size_t gene1, gene2;
    do
    {
        gene1 = randNumber(0, map.size() - 1);
        gene2 = randNumber(0, map.size() - 1);
    } while (gene1 == gene2);
    
  
    const int temp = _genes[gene1];
    _genes[gene1] = _genes[gene2];
    _genes[gene2] = temp;
}

float Individual::getFitness() const {
    return _fitness;
}

Individual Individual::swapGenesOnePoint(const Individual& parent) {
    size_t half_pop = map.size();    
    if (half_pop % 2 != 0)
        half_pop++;
    half_pop = int(half_pop / 2);

    // vector<size_t> parent_genes = parent.getGenesValues();
    vector<size_t> new_genes1;
    vector<size_t> new_genes2;

    for (size_t i = 0; i < half_pop; i++) {
        new_genes1.push_back(_genes[i]);
        // new_genes2.push_back(parent_genes[i]);
    }

    for (size_t i = half_pop; i < map.size(); i++) {
        // new_genes1.push_back(parent_genes[i]);
        new_genes2.push_back(_genes[i]);
    }

    for (size_t i = 0; i < half_pop; i++) {
        new_genes1.push_back(new_genes2[i]);
        // new_genes2.push_back(parent_genes[i]);
    }


    for (size_t i = 0; i < half_pop; i++) {
        new_genes2.push_back(new_genes1[i]);
        // new_genes2.push_back(parent_genes[i]);
    }

    Individual new_indv1(new_genes1);
    new_indv1.calcFitness();
    Individual new_indv2(new_genes2);
    new_indv2.calcFitness();

    // cout << "new gen1: " << new_indv1.getGenes() << endl;
    // cout << "new gen2: " << new_indv2.getGenes() << endl << endl;

    if ( new_indv1.getFitness() < new_indv2.getFitness())
        return new_indv1;    
    else         
        return new_indv2;    
}

vector<size_t> Individual::getGenesValues() const {
    return _genes;
}

string Individual::getGenes() const {
    string str = "";
    for (size_t i = 0; i < map.size(); i++) {
        str += std::to_string(_genes[i]) + " ";
    }
    return str;
}

inline double Individual::calcDistance( double x1, double y1, double x2, double y2) {
    return sqrt( pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0) );
}

inline double Individual::calcRouteDistance() {
    //calcula la suma de todas las distancia entre las rutas sucesivas
    double total_length = 0.0;

    for (size_t i = 0; i < map.size() - 1; i++)
        total_length += calcDistance(
            map[_genes[i]].x,
            map[_genes[i]].y,
            map[_genes[i + 1]].x,
            map[_genes[i + 1]].y );
    // agregamos la distancia entre el primero y el utlimo
    total_length += calcDistance(
            map[_genes[map.size() - 1]].x,
            map[_genes[map.size() - 1]].y,
            map[_genes[0]].x,
            map[_genes[0]].y );

    return total_length;
}


/**************** MAIN ********************/

int main () {

    srand(time(NULL));

    //Coordenadas de prueba
    // int pos_x[] = {0,2,3,5,6,7,4,11,9,2};
    // int pos_y[] = {1,5,1,7,9,3,6,4,0,8};
    //Inicializamos nuestro mapa de coordenadas
    readFile();
    
    // creamos la poblacion inicial de manera aleatoria
    vector<Individual> population;
	for(int i = 0;i<POP_SIZE;i++) {
		population.push_back(Individual());
        population[i].calcFitness();
	}
    
  
    // for(int i = 0;i<POP_SIZE;i++) {
    //     cout << "["<<i<<"]: "<< population[i].getGenes()<< " F: " << population[i].getFitness() <<endl;
	// }
    // cout << "10% " << (10*POP_SIZE)/100 << endl;
    // cout << "90% " << (80*POP_SIZE)/100 << endl;
    int generation = 0;
    float best_fitness;
    // int count = 0;
    while(generation <= 10000) {
        burbleSort(population);
        if (generation == 0)
            float best_fitness = population[0].getFitness();

        if (generation % 1000 == 0) {
            cout<< "Generation: " << generation << endl;
            cout<< "Routes: "<< population[0].getGenes() <<"\t";
            cout<< "Fitness: "<< population[0].getFitness() << endl << endl;
            // if (best_fitness == population[0].getFitness())
            //     count++;
            // if (count == 2 && best_fitness != population[0].getFitness()) {
            //     best_fitness = population[0].getFitness();
            //     count = 0;
            // }
        }

        vector<Individual> new_generation;

        int gen_size = (10*POP_SIZE)/100;
        if (gen_size == 0) gen_size++;
        for (int i = 0; i < gen_size; i++)
            new_generation.push_back(population[i]);

        // gen_size = (90*POP_SIZE)/100;
        for(int i = gen_size;i<POP_SIZE;i+=2) {
            // Individual winner1, winner2;
            // if (count == 2) {
            //     winner1 = tournamentSelection(population, 5, int((20*POP_SIZE)/100));
            //     winner2 = tournamentSelection(population, 5, int((15*POP_SIZE)/100));
            // }
            // else {
            //     winner1 = tournamentSelection(population, 5, int((85*POP_SIZE)/100));
            //     winner2 = tournamentSelection(population, 5, int((75*POP_SIZE)/100));
            // }
            
            Individual winner1 = tournamentSelection(population, 5, int((85*POP_SIZE)/100));
            Individual winner2 = tournamentSelection(population, 5, int((85*POP_SIZE)/100));

            orderCrossover_OX(winner1, winner2);

            winner1.swapMutation();
            winner2.swapMutation();

            winner1.calcFitness();
            winner2.calcFitness();

            new_generation.push_back(winner1);
            new_generation.push_back(winner2);
        }
        population = new_generation;
        generation++;
    }

    return 0;
}