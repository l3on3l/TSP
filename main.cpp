#include <iostream>
using std::cout;
using std::endl;
using std::string;

#include <vector>
using std::vector;

#include <cmath>
using std::sqrt;
using std::pow;

#include <numeric>
#include <algorithm>

#include "main.hpp"

/**************** UTILS ********************/
int randNumber (int start, int end) {
	int range = (end-start)+1;
	int random_int = start+(rand()%range);
	return random_int;
}

bool isGeneInSection(
    vector<size_t> _genes,
    size_t gene,
    size_t iSectionStart,
    size_t iSectionEnd)
{
    for (size_t iGene = iSectionStart; iGene <= iSectionEnd; iGene++)
        if (gene == _genes[iGene])
            return true;

    return false;
}


 vector<size_t> orderCrossover_OX(
    const Individual& dna1,
    const Individual& dna2)
{
    size_t iSecStart = 0;
    size_t iSecEnd = 0;
    
    vector<size_t> _genes;

    while (iSecStart >= iSecEnd)
    {
        iSecStart = randNumber(0, CITY_SIZE - 1);
        iSecEnd   = randNumber(0, CITY_SIZE - 1);
    }
    
    // Copy the section of dna1 to this dna.
    vector<size_t> dna1_genes = dna1.getGenesValues();
    for (size_t iGene = iSecStart; iGene <= iSecEnd; iGene++)
        _genes[iGene] = dna1_genes[iGene];

    const size_t sectionSize = iSecEnd - iSecStart;

    vector<size_t> dna2_genes = dna2.getGenesValues();
    vector<size_t> dnaDifference;
    dnaDifference.reserve(dna2_genes.size() - sectionSize);

    if (iSecEnd + 1 <= dna2_genes.size() - 1)
        for (size_t iGene = iSecEnd + 1; iGene < dna2_genes.size(); iGene++)
            if (!isGeneInSection(_genes, dna2_genes[iGene], iSecStart, iSecEnd))
                dnaDifference.push_back(dna2_genes[iGene]);

    for (size_t iGene = 0; iGene <= iSecEnd; iGene++)
        if (!isGeneInSection(_genes, dna2_genes[iGene], iSecStart, iSecEnd))
            dnaDifference.push_back(dna2_genes[iGene]);
    
    size_t i = 0;

    if (iSecEnd + 1 <= dna2_genes.size() - 1)
        i = iSecEnd + 1;

    for (size_t iGene = 0; iGene < dnaDifference.size(); iGene++)
    {
        _genes[i] = dnaDifference[iGene];
        i++;
        if (i > _genes.size() - 1)
            i = 0;
    }

    return _genes;
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

int getBestFitnessIndex(const  vector<Individual>& population) {
    int best_index = 0;
    for (size_t i = 1; i < population.size()-1; i++) {
        if (population[best_index].getFitness() > population[i].getFitness())
            best_index = i;
    }
    return best_index;
}

Individual tournamentSelection(const vector<Individual>& population, int k) {
    if (k > POP_SIZE) 
        k = POP_SIZE;

    int select[k];
    for (int i = 0; i < k; i++)
        select[i] = randNumber(0, POP_SIZE - 1);

    int bestFitness = select[0];
    for (int i = 1; i < k-1; i++) {
        //el menor fitness es el mejor
        if(population[bestFitness].getFitness() > population[i].getFitness())
            bestFitness = i;
    }
    
    return population[bestFitness];
}

void singlePointCrossovers(const Individual& parent1, const Individual& parent2) {

}

// Individual pmxCrossovers(const Individual& parent1, const Individual& parent2) {
//     int p1 = randNumber(0, CITY_SIZE - 1);
//     int p2 = randNumber(0, CITY_SIZE - 1);
//     int pm = (p1 + p2) / 2;
// }

/**************** CLASS INDIVIDUAL ********************/
Individual::Individual() {
    _genes.resize(CITY_SIZE);
    std::iota(_genes.begin(), _genes.end(), 0);
    initRandGenes();
}

Individual::Individual(const vector<size_t>& genes) {
    _genes.resize(CITY_SIZE);
    for (size_t i = 0; i < _genes.size(); i++)
        _genes[i] = genes[i];
}

Individual::~Individual() { }

void Individual::initRandGenes() {
    for (int i = 0; i < CITY_SIZE; i++)
    {
        size_t rand_gene = randNumber(0, CITY_SIZE - 1);        
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

    int gene1 = randNumber(0, CITY_SIZE-1);
    int gene2 = randNumber(0, CITY_SIZE-1);
  
    const int temp = _genes[gene1];
    _genes[gene1] = _genes[gene2];
    _genes[gene2] = temp;
}

float Individual::getFitness() const {
    return _fitness;
}

Individual Individual::swapGenesOnePoint(const Individual& parent) {
    int half_pop = CITY_SIZE;    
    if (half_pop % 2 != 0)
        half_pop++;
    half_pop = int(half_pop / 2);

    // vector<size_t> parent_genes = parent.getGenesValues();
    vector<size_t> new_genes1;
    vector<size_t> new_genes2;

    for (int i = 0; i < half_pop; i++) {
        new_genes1.push_back(_genes[i]);
        // new_genes2.push_back(parent_genes[i]);
    }

    for (int i = half_pop; i < CITY_SIZE; i++) {
        // new_genes1.push_back(parent_genes[i]);
        new_genes2.push_back(_genes[i]);
    }

    for (int i = 0; i < half_pop; i++) {
        new_genes1.push_back(new_genes2[i]);
        // new_genes2.push_back(parent_genes[i]);
    }


    for (int i = 0; i < half_pop; i++) {
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
    for (int i = 0; i < CITY_SIZE; i++) {
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

    for (int i = 0; i < CITY_SIZE- 1; i++)
        total_length += calcDistance(
            map[_genes[i]].x,
            map[_genes[i]].y,
            map[_genes[i + 1]].x,
            map[_genes[i + 1]].y );
    // agregamos la distancia entre el primero y el utlimo
    total_length += calcDistance(
            map[_genes[CITY_SIZE - 1]].x,
            map[_genes[CITY_SIZE - 1]].y,
            map[_genes[0]].x,
            map[_genes[0]].y );

    return total_length;
}


/**************** MAIN ********************/

int main () {

    srand(time(NULL));

    //Coordenadas de prueba
    int pos_x[] = {0,2,3,5,6,7,4,11,9,2};
    int pos_y[] = {1,5,1,7,9,3,6,4,0,8};
    //Inicializamos nuestro mapa de coordenadas
    for (int i = 0; i < CITY_SIZE; i++) {
        map[i].x = pos_x[i];
        map[i].y = pos_y[i];
    }

    // creamos la poblacion inicial de manera aleatoria
    vector<Individual> population;
	for(int i = 0;i<POP_SIZE;i++) {
		population.push_back(Individual());
        population[i].calcFitness();
	}
    
  
    for(int i = 0;i<POP_SIZE;i++) {
        cout << "["<<i<<"]: "<< population[i].getGenes()<< " F: " << population[i].getFitness() <<endl;
	}

    int generation = 0;
    while(generation < 50) {
        burbleSort(population);
        cout<< "Generation: " << generation + 1 << endl;
        cout<< "Routes: "<< population[0].getGenes() <<"\t";
        cout<< "Fitness: "<< population[0].getFitness() << endl << endl;

        vector<Individual> new_generation;

        int gen_size = (10*POP_SIZE)/100;
        if (gen_size == 0) gen_size++;
        for (int i = 0; i < gen_size; i++)
            new_generation.push_back(population[i]);

        gen_size = (90*POP_SIZE)/100;
        for(int i = 0;i<gen_size;i++) {
            // int i_parent = randNumber(0, gen_size - 1);
            // if (i_parent == i) {
            //     if(i != gen_size - 1)
            //         i_parent++;
            //     else
            //         i_parent--;            
            // }
            Individual winner = tournamentSelection(population, 2);
            winner.calcFitness();
            // Individual new_individual = population[i].swapGenesOnePoint(population[i_parent]);
            Individual new_individual = population[i].swapGenesOnePoint(winner);
            // Individual new_individual(orderCrossover_OX(winner, population[i]));
            new_individual.calcFitness();
            new_individual.swapMutation();
            new_generation.push_back(new_individual);
        }
        population = new_generation;
        generation++;
    }

    return 0;
}