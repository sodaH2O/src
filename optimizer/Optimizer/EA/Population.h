#ifndef __POPULATION_H__
#define __POPULATION_H__

#include <map>
#include <vector>
#include <utility>
#include <queue>
#include <set>
#include <cmath>
#include <sstream>
#include <fstream>

#include "boost/smart_ptr.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

#include "extlib/wfgHypervolume/hypervolume.h"


/**
 *  \class Population
 *
 *  Managing a population of individuals. We maintain two sets: a set of all
 *  (evaluated) individuals in the population and a set of new potential
 *  individuals (the selector decides which individuals join the population),
 *  called 'stagingArea'.
 *  Most operations work on the 'stagingArea', population is kept for
 *  visualization purposes.
 */
template< class Individual_t >
class Population {

public:
    Population() {
        last_identity = 0;
        hypervolume_  = 0.0;
    }

    ~Population() {
        clean_population();
    }

    typedef typename Individual_t::genes_t          genes_t;
    typedef boost::shared_ptr<Individual_t>         individual;
    typedef std::pair< unsigned int, individual >   ind_t;

    /// population iterator type
    typedef typename std::map<unsigned int, individual>::iterator indItr_t;

    /**
     *  Adds an individual to the population
     *  @param ind an individual that will be added to the population
     *  @return individual id if successful, -1 otherwise
     */
    unsigned int add_individual(individual ind) {

        unsigned int id = getFreeID();
        stagingArea.insert(ind_t(id, ind));
        ind->id = id;
        //std::cout << "+++ staging    id = " << id << "\xd";
        return id;
    }


    void remove_individual(individual ind) {

        indItr_t it = stagingArea.begin();
        while(it != stagingArea.end()) {
            if(it->second == ind) {
                if(it->first == last_identity-1)
                    last_identity--;
                else
                    freeids.push(it->first);

                //std::cout << "--- removing   id = " << it->first << "\xd";
                stagingArea.erase(it);
                break;
            }
            it++;
        }
    }


    /**
     *  Get an individual of the current population with a specific ID
     *  @param identity an ID of the individual in the population
     *  @return the individual with the specified ID in the population, -1 if
     *  none found
     */
    individual get_individual(int identity) {

        indItr_t it;
        if(identity == -1)
            it = individuals.begin();
        else {
            it = individuals.find(identity);
            if( it == individuals.end() )
                return individual();
        }

        return it->second;
    }


    individual get_staging(int identity) {

        indItr_t it;
        if(identity == -1)
            it = stagingArea.begin();
        else {
            it = stagingArea.find(identity);
            if( it == stagingArea.end() )
                return individual();
        }

        return it->second;
    }


    void commit_individuals(std::set<unsigned int> ids) {

        foreach(unsigned int id, ids) {
            //std::cout << "--+ committing id = " << id << "\xd";
            individual ind = get_staging(id);
            individuals.insert(ind_t(id, ind));
            stagingArea.erase(id);
        }
    }

    /**
     *  Remove all non-surviving individuals from the population and put ID's
     *  back in pool of free ID's.
     *  @param surviors to keep for next generation
     */
    void keepSurvivors(std::set<unsigned int> survivors) {

        indItr_t it = individuals.begin();
        while(it != individuals.end()) {
            if( survivors.count(it->first) == 0 ) {
                if(it->first == last_identity-1)
                    last_identity--;
                else
                    freeids.push(it->first);

                individuals.erase(it++);
            } else
                it++;
        }
    }


    /// check if a gene set is already represented in the population
    //XXX: currently O(n): add a fast look-up table?
    bool isRepresentedInPopulation(genes_t ind_genes) {

        foreach(ind_t ind, individuals) {
            if( ind_genes == ind.second->genes )
                return true;
        }

        return false;
    }


    double computeHypervolume(size_t island_id) {

        std::ofstream file;
        std::ostringstream filename;
        filename << "hypervol.dat_" << island_id;
        file.open(filename.str().c_str(), std::ios::out);

        file << "#" << std::endl;

        indItr_t it;
        for(it = individuals.begin(); it != individuals.end(); it++) {

            individual temp = it->second;
            for(size_t i=0; i<temp->objectives.size(); i++)
                file << temp->objectives[i] << " ";
            if (temp->objectives.size() > 0)
                file << std::endl;
        }

        file << "#" << std::endl;

        file.flush();
        file.close();

        hypervolume_ = Hypervolume::FromFile(filename.str());
        return hypervolume_;
    }


    void commit_individuals() {
        individuals.insert(stagingArea.begin(), stagingArea.end());
        stagingArea.clear();
    }


    /// iterator begin on staging area
    indItr_t stagingBegin() { return stagingArea.begin(); }
    /// iterator end on staging area
    indItr_t stagingEnd()   { return stagingArea.end(); }


    /**
     *  Size of population
     *  @return total number of individuals in population
     */
    unsigned int size() const { return individuals.size(); }

    /// iterator begin on population container
    indItr_t begin() { return individuals.begin(); }
    /// iterator end on population container
    indItr_t end()   { return individuals.end(); }

    /// Remove (and cleanup) all individuals in the population
    void clean_population() {

        stagingArea.clear();
        individuals.clear();
    }

private:

    /// population container holding all individuals
    std::map<unsigned int, individual > individuals;

    /// staging area for individuals probably joining population
    std::map<unsigned int, individual > stagingArea;

    /// queue to handle free individual ID's
    std::queue<unsigned int> freeids;

    /// last used (= next free) ID
    unsigned int last_identity;

    double hypervolume_;

    /**
     *  Manages free individual ID's
     *  @return lowest free ID
     */
    unsigned int getFreeID() {

        unsigned int id = 0;
        if(freeids.empty()) {
            id = last_identity;
            last_identity++;
        } else {
            id = freeids.front();
            freeids.pop();
        }

        return id;
    }
};

#endif