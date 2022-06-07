#include<vector>
using std::vector;

#include "individual.hpp"
class GA {
    private:
        vector<Individual> _population;
        double _total_fitness;
        int _generation;
    public:
        //propios
        GA();
        ~GA();
        //getters & setters
        double getTotalFitness() const;
        void setPopulation();
        //metodos
        void calcTotalFitness();
        void proportionalFractionPop();
        void updateGeneration();
        void fitnessProportionateSelection();
        Individual crossovers(const Individual& parent1, const Individual& parent2);
};