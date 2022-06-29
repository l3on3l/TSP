#include <iostream>
using std::string;

#include <vector>
using std::vector;
class Individual {    
    private:
        vector<size_t> _genes;
        double _fitness;
        double _proportional_faction;
    public:
        //propios
        Individual();
        Individual(const vector<size_t>& genes);
        ~Individual();
        //getters & setters
        float getFitness() const;
        string getGenes() const;
        vector<size_t> getGenesValues() const;
        Individual swapGenesOnePoint(const Individual& parent);
        void swapMutation();
        void setProportionalFraction(const double value);
        //metodos
        void calcFitness();
        Individual mate(const Individual& parent);        
        void mutateGenes(double probability); //falta
        void initRandGenes();
        bool operator<(Individual &p);
        //friends
        friend void orderCrossover_OX(Individual& parent1, Individual& parent2);
    private:
        inline double calcDistance(double,double,double,double);
        inline double calcRouteDistance();
};