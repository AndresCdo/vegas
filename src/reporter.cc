#include "../include/reporter.h"
#include "../include/exception.h"

#include <array>

Reporter::Reporter()
{
    initializeHandles();
}

void Reporter::initializeHandles()
{
    this->file = -1;
    this->dataspace_id_energy = -1;
    this->memspace_id_ = -1;
    this->status = 0;
    this->energies_dset = -1;
    this->temps_dset = -1;
    this->fields_dset = -1;
    this->position_dset = -1;
    this->types_dset = -1;
    this->finalstates_dset = -1;
    this->memspace_id_finalstates = -1;
    this->dataspace_id_finalstates = -1;
    this->closed_ = false;
}

hid_t Reporter::createCompressionProperties()
{
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0) {
        if (this->file >= 0) {
            H5Fclose(this->file);
            this->file = -1;
        }
        throw vegas::HDF5Exception("Failed to create HDF5 property list");
    }
    return dcpl;
}

void Reporter::createMagnetizationDatasets(Lattice& lattice, Index numTemps, Index mcs, hid_t dcpl)
{
    Index num_types = lattice.getMapTypeIndexes().size();
    
    this->mags_dset_x_ = std::vector<hid_t>(num_types + 1, -1);
    this->dataspace_id_mag_x_ = std::vector<hid_t>(num_types + 1, -1);
    this->mags_dset_y_ = std::vector<hid_t>(num_types + 1, -1);
    this->dataspace_id_mag_y_ = std::vector<hid_t>(num_types + 1, -1);
    this->mags_dset_z_ = std::vector<hid_t>(num_types + 1, -1);
    this->dataspace_id_mag_z_ = std::vector<hid_t>(num_types + 1, -1);

    hsize_t dims[2] = {numTemps, mcs};
    hid_t space_mag = H5Screate_simple(2, dims, NULL);
    if (space_mag < 0) {
        throw vegas::HDF5Exception("Failed to create HDF5 dataspace for magnetization");
    }

    for (auto& type : lattice.getMapTypeIndexes())
    {
        this->mags_dset_x_.at(type.second) = H5Dcreate(file, (type.first + "_x").c_str(),
                    H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        this->mags_dset_y_.at(type.second) = H5Dcreate(file, (type.first + "_y").c_str(),
                    H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        this->mags_dset_z_.at(type.second) = H5Dcreate(file, (type.first + "_z").c_str(),
                    H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    }

    this->mags_dset_x_.at(num_types) = H5Dcreate(file, "magnetization_x",
                H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    this->mags_dset_y_.at(num_types) = H5Dcreate(file, "magnetization_y",
                H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    this->mags_dset_z_.at(num_types) = H5Dcreate(file, "magnetization_z",
                H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    this->energies_dset = H5Dcreate(file, "energy",
                H5T_IEEE_F64LE, space_mag, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    this->status = H5Sclose(space_mag);
}

void Reporter::createAuxiliaryDatasets(Lattice& lattice, Index numTemps)
{
    hsize_t dims_temps[1] = {numTemps};
    hid_t space_temp = H5Screate_simple(1, dims_temps, NULL);
    if (space_temp >= 0) {
        this->temps_dset = H5Dcreate(file, "temperature",
                    H5T_IEEE_F64LE, space_temp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        this->status = H5Sclose(space_temp);
    }

    hsize_t dims_fields[1] = {numTemps};
    hid_t space_field = H5Screate_simple(1, dims_fields, NULL);
    if (space_field >= 0) {
        this->fields_dset = H5Dcreate(file, "field",
                    H5T_IEEE_F64LE, space_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        this->status = H5Sclose(space_field);
    }

    hsize_t dims_pos[2] = {lattice.getAtoms().size(), 3};
    hid_t space_pos = H5Screate_simple(2, dims_pos, NULL);
    if (space_pos >= 0) {
        this->position_dset = H5Dcreate(file, "positions",
                    H5T_IEEE_F64LE, space_pos, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        this->status = H5Sclose(space_pos);
    }

    hid_t filetype = H5Tcopy(H5T_C_S1);
    this->status = H5Tset_size(filetype, H5T_VARIABLE);
    
    hsize_t dims_types[1] = {lattice.getAtoms().size()};
    hid_t space_types = H5Screate_simple(1, dims_types, NULL);
    if (space_types >= 0) {
        this->types_dset = H5Dcreate(file, "types",
                    filetype, space_types, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        this->status = H5Sclose(space_types);
    }
    
    this->status = H5Tclose(filetype);

    hsize_t dims_finalstates[3] = {numTemps, lattice.getAtoms().size(), 3};
    hid_t space_final = H5Screate_simple(3, dims_finalstates, NULL);
    if (space_final >= 0) {
        this->finalstates_dset = H5Dcreate(file, "finalstates",
                    H5T_IEEE_F64LE, space_final, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        this->status = H5Sclose(space_final);
    }
}

void Reporter::writeStaticData(Lattice& lattice, const std::vector<Real>& temps, const std::vector<Real>& fields)
{
    std::vector<std::array<double, 3>> positions(lattice.getAtoms().size());
    std::vector<const char*> types(lattice.getAtoms().size());
    
    Index i = 0;
    for (auto& atom : lattice.getAtoms())
    {
        std::copy(std::begin(atom.getPosition()), std::end(atom.getPosition()), positions[i].begin());
        types[i] = atom.getType().c_str();
        i++;
    }

    hid_t memtype_write = H5Tcopy(H5T_C_S1);
    this->status = H5Tset_size(memtype_write, H5T_VARIABLE);
    
    this->status = H5Dwrite(types_dset, memtype_write, H5S_ALL, H5S_ALL, H5P_DEFAULT, types.data());
    this->status = H5Dwrite(position_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions.data()->data());
    this->status = H5Dwrite(temps_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temps.data());
    this->status = H5Dwrite(fields_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fields.data());
    
    if (memtype_write >= 0) {
        this->status = H5Tclose(memtype_write);
    }
}

void Reporter::setupHyperslabs(Index mcs, Index numTemps, Index numAtoms)
{
    (void)numTemps;
    (void)numAtoms;
    
    this->dims_select_[0] = mcs;
    this->memspace_id_ = H5Screate_simple(1, this->dims_select_, NULL);

    this->dataspace_id_energy = H5Dget_space(this->energies_dset);

    this->start_[1] = 0;
    this->count_[0] = 1;
    this->count_[1] = mcs;
    this->stride_[0] = 1;
    this->stride_[1] = 1;
    this->block_[0] = 1;
    this->block_[1] = 1;

    Index num_types = static_cast<Index>(mags_dset_x_.size() - 1);
    for (Index idx = 0; idx <= num_types; ++idx)
    {
        if (mags_dset_x_.at(idx) >= 0) {
            this->dataspace_id_mag_x_.at(idx) = H5Dget_space(mags_dset_x_.at(idx));
        }
        if (mags_dset_y_.at(idx) >= 0) {
            this->dataspace_id_mag_y_.at(idx) = H5Dget_space(mags_dset_y_.at(idx));
        }
        if (mags_dset_z_.at(idx) >= 0) {
            this->dataspace_id_mag_z_.at(idx) = H5Dget_space(mags_dset_z_.at(idx));
        }
    }

    this->dims_select_finalstates[0] = 3;
    this->memspace_id_finalstates = H5Screate_simple(1, this->dims_select_finalstates, NULL);
    this->dataspace_id_finalstates = H5Dget_space(this->finalstates_dset);

    this->start_finalstates[2] = 0;
    this->count_finalstates[0] = 1;
    this->count_finalstates[1] = 1;
    this->count_finalstates[2] = 3;
    this->stride_finalstates[0] = 1;
    this->stride_finalstates[1] = 1;
    this->stride_finalstates[2] = 1;
    this->block_finalstates[0] = 1;
    this->block_finalstates[1] = 1;
    this->block_finalstates[2] = 1;
}

void Reporter::writeAttributes(Index mcs, Index seed, Real kb)
{
    hid_t aid2 = H5Screate(H5S_SCALAR);
    if (aid2 >= 0) {
        hid_t attr_mcs = H5Acreate(file, "mcs", H5T_NATIVE_INT, aid2, H5P_DEFAULT, H5P_DEFAULT);
        if (attr_mcs >= 0) {
            this->status = H5Awrite(attr_mcs, H5T_NATIVE_INT, &mcs);
            this->status = H5Aclose(attr_mcs);
        }
        
        hid_t attr_seed = H5Acreate(file, "seed", H5T_NATIVE_INT, aid2, H5P_DEFAULT, H5P_DEFAULT);
        if (attr_seed >= 0) {
            this->status = H5Awrite(attr_seed, H5T_NATIVE_INT, &seed);
            this->status = H5Aclose(attr_seed);
        }
        
        hid_t attr_kb = H5Acreate(file, "kb", H5T_NATIVE_DOUBLE, aid2, H5P_DEFAULT, H5P_DEFAULT);
        if (attr_kb >= 0) {
            this->status = H5Awrite(attr_kb, H5T_NATIVE_DOUBLE, &kb);
            this->status = H5Aclose(attr_kb);
        }
        
        this->status = H5Sclose(aid2);
    }
}

Reporter::Reporter(std::string filename,
             std::vector<Array> magnetizationTypes,
             Lattice& lattice,
             const std::vector<Real>& temps,
             const std::vector<Real>& fields,
             Index mcs,
             Index seed,
             Real kb)
{
    (void)magnetizationTypes;  // Unused parameter
    
    initializeHandles();

    this->file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (this->file < 0) {
        this->file = -1;
        throw vegas::HDF5Exception("Failed to create HDF5 file: " + filename);
    }

    hid_t dcpl = createCompressionProperties();

    hsize_t CHUNK[2] = {1, Index(mcs / AMOUNTCHUNKS)};
    this->status = H5Pset_deflate(dcpl, 1);
    this->status = H5Pset_chunk(dcpl, 2, CHUNK);

    createMagnetizationDatasets(lattice, temps.size(), mcs, dcpl);
    createAuxiliaryDatasets(lattice, temps.size());
    
    this->status = H5Pclose(dcpl);

    writeStaticData(lattice, temps, fields);
    setupHyperslabs(mcs, temps.size(), lattice.getAtoms().size());
    writeAttributes(mcs, seed, kb);
}


void Reporter::partial_report(
    const std::vector<Real>& enes,
    const std::vector< std::vector<Real> >& histMag_x,
    const std::vector< std::vector<Real> >& histMag_y,
    const std::vector< std::vector<Real> >& histMag_z,
    Lattice& lattice, Index index)
{
    this->start_[0] = index;

    this->status = H5Sselect_hyperslab(this->dataspace_id_energy, H5S_SELECT_SET, this->start_,
                                  this->stride_, this->count_, this->block_);
    this->status = H5Dwrite(this->energies_dset, H5T_NATIVE_DOUBLE, this->memspace_id_,
                       this->dataspace_id_energy, H5P_DEFAULT, enes.data());

    Index i = 0;
    for (auto& val : this->mags_dset_x_)
    {
        this->status = H5Sselect_hyperslab(this->dataspace_id_mag_x_.at(i), H5S_SELECT_SET, this->start_,
                                      this->stride_, this->count_, this->block_);
        this->status = H5Dwrite(val, H5T_NATIVE_DOUBLE, this->memspace_id_,
                                   this->dataspace_id_mag_x_.at(i), H5P_DEFAULT, histMag_x.at(i).data());

        this->status = H5Sselect_hyperslab(this->dataspace_id_mag_y_.at(i), H5S_SELECT_SET, this->start_,
                                      this->stride_, this->count_, this->block_);
        this->status = H5Dwrite(mags_dset_y_.at(i), H5T_NATIVE_DOUBLE, this->memspace_id_,
                                   this->dataspace_id_mag_y_.at(i), H5P_DEFAULT, histMag_y.at(i).data());

        this->status = H5Sselect_hyperslab(this->dataspace_id_mag_z_.at(i), H5S_SELECT_SET, this->start_,
                                      this->stride_, this->count_, this->block_);
        this->status = H5Dwrite(mags_dset_z_.at(i), H5T_NATIVE_DOUBLE, this->memspace_id_,
                                   this->dataspace_id_mag_z_.at(i), H5P_DEFAULT, histMag_z.at(i).data());
        i++;
    }

    this->start_finalstates[0] = index;
    i = 0;
    for (auto& atom : lattice.getAtoms())
    {
        this->start_finalstates[1] = i;
        std::vector<double> spin;
        spin.assign(std::begin(atom.getSpin()), std::end(atom.getSpin()));

        this->status = H5Sselect_hyperslab(this->dataspace_id_finalstates, H5S_SELECT_SET, this->start_finalstates,
                                             this->stride_finalstates, this->count_finalstates, this->block_finalstates);
        this->status = H5Dwrite(this->finalstates_dset, H5T_NATIVE_DOUBLE, this->memspace_id_finalstates,
                                   this->dataspace_id_finalstates, H5P_DEFAULT, spin.data());

        i++;
    }
}

void Reporter::close()
{
    if (this->closed_) {
        return;
    }

    for (size_t i = 0; i < dataspace_id_mag_x_.size(); ++i)
    {
        if (dataspace_id_mag_x_[i] >= 0) {
            this->status = H5Sclose(dataspace_id_mag_x_[i]);
            dataspace_id_mag_x_[i] = -1;
        }
        if (dataspace_id_mag_y_[i] >= 0) {
            this->status = H5Sclose(dataspace_id_mag_y_[i]);
            dataspace_id_mag_y_[i] = -1;
        }
        if (dataspace_id_mag_z_[i] >= 0) {
            this->status = H5Sclose(dataspace_id_mag_z_[i]);
            dataspace_id_mag_z_[i] = -1;
        }
    }

    if (dataspace_id_energy >= 0) {
        this->status = H5Sclose(dataspace_id_energy);
        dataspace_id_energy = -1;
    }
    if (memspace_id_ >= 0) {
        this->status = H5Sclose(memspace_id_);
        memspace_id_ = -1;
    }
    if (memspace_id_finalstates >= 0) {
        this->status = H5Sclose(memspace_id_finalstates);
        memspace_id_finalstates = -1;
    }
    if (dataspace_id_finalstates >= 0) {
        this->status = H5Sclose(dataspace_id_finalstates);
        dataspace_id_finalstates = -1;
    }

    Index i = 0;
    for (auto& val : this->mags_dset_x_)
    {
        if (val >= 0) {
            this->status = H5Dclose(val);
            val = -1;
        }
        if (mags_dset_y_.at(i) >= 0) {
            this->status = H5Dclose(mags_dset_y_.at(i));
            mags_dset_y_.at(i) = -1;
        }
        if (mags_dset_z_.at(i) >= 0) {
            this->status = H5Dclose(mags_dset_z_.at(i));
            mags_dset_z_.at(i) = -1;
        }
        i++;
    }

    if (energies_dset >= 0) {
        this->status = H5Dclose(this->energies_dset);
        energies_dset = -1;
    }
    if (temps_dset >= 0) {
        this->status = H5Dclose(this->temps_dset);
        temps_dset = -1;
    }
    if (fields_dset >= 0) {
        this->status = H5Dclose(this->fields_dset);
        fields_dset = -1;
    }
    if (position_dset >= 0) {
        this->status = H5Dclose(this->position_dset);
        position_dset = -1;
    }
    if (types_dset >= 0) {
        this->status = H5Dclose(this->types_dset);
        types_dset = -1;
    }
    if (finalstates_dset >= 0) {
        this->status = H5Dclose(this->finalstates_dset);
        finalstates_dset = -1;
    }

    if (file >= 0) {
        this->status = H5Fclose(this->file);
        file = -1;
    }

    this->closed_ = true;
}

Reporter::~Reporter()
{
    close();
}
