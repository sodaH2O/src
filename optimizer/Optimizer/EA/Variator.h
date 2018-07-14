#ifndef __VARIATOR_H__
#define __VARIATOR_H__

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Optimizer/EA/Population.h"
#include "Optimizer/Optimizer.h"

template<
      class ind_t
    , template <class> class CrossoverOperator
    , template <class> class MutationOperator
>
class Variator : public CrossoverOperator<ind_t>,
                 public MutationOperator<ind_t>
{

public:

    Variator(size_t sizeInitial,
             std::vector<std::string> dNames,
             Optimizer::bounds_t dVarBounds, Expressions::Named_t constraints,
             CmdArguments_t args)
        : sizeInitial_m(sizeInitial)
        , dNames_m(dNames)
        , dVarBounds_m(dVarBounds)
    {
        // add constraints, if only design variables are needed for evaluation
        for(auto constraint : constraints) {
            bool allDesignVariables = true;
            std::set<std::string> req_vars = constraint.second->getReqVars();
            if (req_vars.empty()) allDesignVariables = false;
            for (std::string req_var : req_vars) {
                // check if it is a design variable
                if (std::find(dNames_m.begin(),dNames_m.end(),req_var) == dNames_m.end()) {
                    allDesignVariables = false;
                    break;
                }
            }
            if (allDesignVariables == true)
                constraints_m.insert(constraint);
        }

        //FIXME: pass population as arg to variator
        //boost::shared_ptr< Population<ind_t> >
        population_m.reset(new Population<ind_t>());

        mutationProbability_m =
            args->getArg<double>("mutation-probability", 0.5);

        recombinationProbability_m =
            args->getArg<double>("recombination-probability", 0.5);

        args_ = args;
    }

    ~Variator() {
    }

    //FIXME access population from outside
    boost::shared_ptr< Population<ind_t> > population() {
        return population_m;
    }

    /// create an initial population
    void initial_population() {
        for(size_t i = 0; i < sizeInitial_m; i++)
            new_individual();

    }

    /// set an individual as individual: replace with a new individual
    void infeasible(boost::shared_ptr<ind_t> ind) {
        population_m->remove_individual(ind);
        new_individual();
    }

    /// returns false if all individuals have been evaluated
    bool hasMoreIndividualsToEvaluate() {
        return !individualsToEvaluate_m.empty();
    }

    /// return next individual to evaluate
    boost::shared_ptr<ind_t> popIndividualToEvaluate() {
        unsigned int ind = individualsToEvaluate_m.front();
        individualsToEvaluate_m.pop();
        return population_m->get_staging(ind);
    }

    /** Performs variation (recombination and mutation) on a set of parent
     *  individuals.
     *
     *  @param[in] parents
     */
    void variate(std::vector<unsigned int> parents) {

        // copying all individuals from parents
        for(unsigned int parent : parents) {
            new_individual( population_m->get_individual(parent) );
        }

        // only variate new offspring, individuals in staging area have been
        // variated already
        std::queue<unsigned int> tmp(individualsToEvaluate_m);
        while(!tmp.empty()) {

            // pop first individual
            unsigned int idx = tmp.front(); tmp.pop();
            boost::shared_ptr<ind_t> a = population_m->get_staging(idx);

            // handle special case where we have an uneven number of offspring
            if(tmp.empty()) {
                if (drand(1) <= mutationProbability_m)
                    this->mutate(a, args_);
                break;
            }

            // and second if any
            idx = tmp.front(); tmp.pop();
            boost::shared_ptr<ind_t> b = population_m->get_staging(idx);

            // create new individuals

            // temporary copy in case not successful
            boost::shared_ptr<ind_t> copyA(new ind_t(a));
            boost::shared_ptr<ind_t> copyB(new ind_t(b));

            int iter = 0;
            while (true) {
                // assign with shared pointer constructor
                *a = copyA;
                *b = copyB;

                // do recombination
                if(drand(1) <= recombinationProbability_m) {
                    this->crossover(a, b, args_);
                }

                // do mutation
                if (drand(1) <= mutationProbability_m) {
                    this->mutate(a, args_);
                    this->mutate(b, args_);
                }

                // check if viable offspring
                bool viableA = a->viable();
                bool viableB = b->viable();
                if (viableA == true && viableB == true) {
                    break;
                }
		std::cout << "Individual not viable, I try again: iter= " << iter << std::endl;
                iter++;
                // if maximum number of tries then create new individual(s)
                if (iter > 100) {
                    if (viableA == false) {
                        infeasible(a);
                    }
                    if (viableB == false) {
                        infeasible(b);
                    }
                    break;
                }
            }
        }
    }


protected:

    /// create a new individual
    void new_individual() {
        boost::shared_ptr<ind_t> ind(new ind_t(dVarBounds_m, dNames_m, constraints_m));
        individualsToEvaluate_m.push( population_m->add_individual(ind) );
    }

    /// copy an individual
    void new_individual(boost::shared_ptr<ind_t> ind) {
        boost::shared_ptr<ind_t> return_ind(new ind_t(ind));
        individualsToEvaluate_m.push(
            population_m->add_individual(return_ind) ) ;
    }

private:

    /// population of individuals
    boost::shared_ptr< Population<ind_t> > population_m;
    /// number of individuals in initial population
    size_t sizeInitial_m;

    /// user specified command line arguments
    CmdArguments_t args_;

    /// keep a queue of individuals that have to be evaluated
    std::queue<unsigned int> individualsToEvaluate_m;

    /// names of the design variables
    std::vector<std::string> dNames_m;
    /// bounds on design variables
    Optimizer::bounds_t dVarBounds_m;
    /// constraints
    Expressions::Named_t constraints_m;

    /// probability of applying the mutation operator
    double mutationProbability_m;
    /// probability of applying the recombination operator
    double recombinationProbability_m;

    /**
     *  Get a random double between [0, range]
     *  @param[in] range of random number
     *  @return random double value between [0, range]
     */
    double drand(double range) {
        return (range * (double) rand() / (RAND_MAX + 1.0));
    }

    /**
     *  Get a random integer between [0, range]
     *  @param[in] range of random number
     *  @return random integer value between [0, range]
     */
    int irand(int range) {
        return (int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
    }

};

#endif
