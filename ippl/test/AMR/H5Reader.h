#ifndef H5READER_H
#define H5READER_H

#include <string>

#include <H5hut.h>

#include "Distribution.h"
#include "helper_functions.h" // container_t

#define READDATA(type, file, name, value) H5PartReadData##type(file, name, value);
#define WRITEDATA(type, file, name, value) H5PartWriteData##type(file, name, value);
#define WRITESTRINGFILEATTRIB(file, name, value) H5WriteFileAttribString(file, name, value);

/*!
 * @file H5Reader.h
 * @details This class is implemented according to
 * src/Structure/H5PartWrapper and
 * src/Structure/H5PartWrapperForPC.\n It is used in
 * Distribution.h to read in a particle distribution
 * from a H5 file (supports only opal flavour: opal-cycl).
 * @author Matthias Frey
 * @date 25. - 26. October 2016
 * @brief Read in a particle distribution from a H5 file
 */

/// Defines functions for reading a distribution from a cyclotron H5 file
class H5Reader {
    
public:
    /*!
     * @param filename specifies path and name of the file
     */
    H5Reader(const std::string& filename);
    
    H5Reader();
    
    /*!
     * Open the file and set the step to read / write
     * @param step to read in / write
     * @param flags (read) H5_O_RDONLY or (write) H5_O_WRONLY
     */
    void open(int step, h5_int32_t flags = H5_O_RDONLY);
    
    /*!
     * Close the file and sets the pointer to NULL
     */
    void close();
    
    /*!
     * Copy the particle distribution to the containers of
     * the Distribution class. (parallel read)
     * It performs a shift in longitduinal direction due to the
     * box extent
     * @param x - coordinate
     * @param px - coordinate
     * @param y - coordinate
     * @param py - coordinate
     * @param z - coordinate
     * @param pz - coordinate
     * @param q is the particle charge
     * @param mass of a particle
     * @param firstParticle to read (core specific)
     * @param lastParticle to read (core specific)
     */
    void read(Distribution::container_t& x,
              Distribution::container_t& px,
              Distribution::container_t& y,
              Distribution::container_t& py,
              Distribution::container_t& z,
              Distribution::container_t& pz,
              Distribution::container_t& q,
              Distribution::container_t& mass,
              size_t firstParticle,
              size_t lastParticle);
    
    
#ifdef IPPL_AMR
    void writeHeader();
    
    /*!
     * Write a particle distribution to a file
     */
    void write(PartBunchAmr< ParticleAmrLayout<double, AMREX_SPACEDIM> >* bunch);
#endif
    
    /*!
     * @returns the number of particles
     */
    h5_ssize_t getNumParticles();
    
    void writeScalarField(const container_t& scalfield,
                          const Array<Geometry>& geom);
    
    void writeVectorField(const container_t& vecfield);
    
private:
    std::string filename_m;     ///< Path and filename
    h5_file_t file_m;           ///< Opened file
};

#endif
