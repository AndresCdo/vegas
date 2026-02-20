#ifndef REPORTER
#define REPORTER

#include "params.h"
#include "lattice.h"
#include <hdf5.h>

#include <string>
#include <vector>
#include <map>

class Reporter
{
public:
    Reporter();
    Reporter(std::string filename,
             std::vector<Array> magnetizationTypes,
             Lattice& lattice,
             const std::vector<Real>& temps,
             const std::vector<Real>& fields,
             Index mcs,
             Index seed,
             Real kb);

    void partial_report(
        const std::vector<Real>& enes,
        const std::vector< std::vector<Real> >& histMag_x,
        const std::vector< std::vector<Real> >& histMag_y,
        const std::vector< std::vector<Real> >& histMag_z,
        Lattice& lattice, Index index);
    void close();
    ~Reporter();

private:
    hid_t       file;
    hid_t       dataspace_id_energy, memspace_id_;
    herr_t      status;

    std::vector<hid_t> mags_dset_x_;
    std::vector<hid_t> dataspace_id_mag_x_;

    std::vector<hid_t> mags_dset_y_;
    std::vector<hid_t> dataspace_id_mag_y_;

    std::vector<hid_t> mags_dset_z_;
    std::vector<hid_t> dataspace_id_mag_z_;

    hid_t energies_dset;
    hid_t temps_dset;
    hid_t fields_dset;

    hid_t position_dset;
    hid_t types_dset;
    hid_t finalstates_dset;

    hsize_t     count_[2];
    hsize_t     start_[2];
    hsize_t     stride_[2];
    hsize_t     block_[2];

    hsize_t     dims_select_[1];

    hid_t       memspace_id_finalstates, dataspace_id_finalstates;
    hsize_t     count_finalstates[3];
    hsize_t     start_finalstates[3];
    hsize_t     stride_finalstates[3];
    hsize_t     block_finalstates[3];

    hsize_t     dims_select_finalstates[1];

    // Private helper methods for constructor initialization
    void initializeHandles();
    hid_t createCompressionProperties();
    void createMagnetizationDatasets(Lattice& lattice, Index numTemps, Index mcs, hid_t dcpl);
    void createAuxiliaryDatasets(Lattice& lattice, Index numTemps);
    void writeStaticData(Lattice& lattice, const std::vector<Real>& temps, const std::vector<Real>& fields);
    void setupHyperslabs(Index mcs, Index numTemps, Index numAtoms);
    void writeAttributes(Index mcs, Index seed, Real kb);
};

#endif
