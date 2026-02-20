// Reporter unit tests for VEGAS simulation package

#include "../include/reporter.h"
#include "../include/lattice.h"
#include "../include/params.h"
#include <iostream>
#include <cassert>
#include <cstdio>
#include <hdf5.h>

void test_reporter_default_constructor()
{
    std::cout << "Testing Reporter default constructor... ";
    
    Reporter reporter;
    
    std::cout << "PASS\n";
}

void test_reporter_constructor_creates_file()
{
    std::cout << "Testing Reporter constructor creates file... ";
    
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    std::vector<Real> temps = {1.0};
    std::vector<Real> fields = {0.0};
    std::vector<Array> magTypes;
    
    Reporter reporter("test_reporter_output.h5", magTypes, lattice, temps, fields, 10, 42, 1.0);
    
    FILE* f = fopen("test_reporter_output.h5", "r");
    assert(f != nullptr);
    fclose(f);
    
    std::remove("test_reporter_output.h5");
    
    std::cout << "PASS\n";
}

void test_reporter_file_structure()
{
    std::cout << "Testing Reporter file structure... ";
    
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    std::vector<Real> temps = {1.0, 2.0};
    std::vector<Real> fields = {0.0, 1.0};
    std::vector<Array> magTypes;
    
    Reporter reporter("test_structure_output.h5", magTypes, lattice, temps, fields, 10, 42, 1.0);
    
    hid_t file = H5Fopen("test_structure_output.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    
    hid_t dset = H5Dopen(file, "magnetization_x", H5P_DEFAULT);
    assert(dset >= 0);
    H5Dclose(dset);
    
    dset = H5Dopen(file, "energy", H5P_DEFAULT);
    assert(dset >= 0);
    H5Dclose(dset);
    
    dset = H5Dopen(file, "temperature", H5P_DEFAULT);
    assert(dset >= 0);
    H5Dclose(dset);
    
    dset = H5Dopen(file, "positions", H5P_DEFAULT);
    assert(dset >= 0);
    H5Dclose(dset);
    
    H5Fclose(file);
    
    std::remove("test_structure_output.h5");
    
    std::cout << "PASS\n";
}

void test_reporter_attributes()
{
    std::cout << "Testing Reporter attributes... ";
    
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    std::vector<Real> temps = {1.0};
    std::vector<Real> fields = {0.0};
    std::vector<Array> magTypes;
    
    Reporter reporter("test_attrs_output.h5", magTypes, lattice, temps, fields, 100, 12345, 1.5);
    
    hid_t file = H5Fopen("test_attrs_output.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    
    hid_t attr = H5Aopen(file, "mcs", H5P_DEFAULT);
    assert(attr >= 0);
    
    int mcs_val;
    H5Aread(attr, H5T_NATIVE_INT, &mcs_val);
    assert(mcs_val == 100);
    H5Aclose(attr);
    
    attr = H5Aopen(file, "seed", H5P_DEFAULT);
    assert(attr >= 0);
    int seed_val;
    H5Aread(attr, H5T_NATIVE_INT, &seed_val);
    assert(seed_val == 12345);
    H5Aclose(attr);
    
    attr = H5Aopen(file, "kb", H5P_DEFAULT);
    assert(attr >= 0);
    double kb_val;
    H5Aread(attr, H5T_NATIVE_DOUBLE, &kb_val);
    assert(fp_equal(kb_val, 1.5));
    H5Aclose(attr);
    
    H5Fclose(file);
    
    std::remove("test_attrs_output.h5");
    
    std::cout << "PASS\n";
}

void test_reporter_partial_report()
{
    std::cout << "Testing Reporter partial report... ";
    
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    std::vector<Real> temps = {1.0};
    std::vector<Real> fields = {0.0};
    std::vector<Array> magTypes;
    
    // Get number of types to create proper history vectors
    Index numTypes = lattice.getMapTypeIndexes().size();
    
    Reporter reporter("test_partial_output.h5", magTypes, lattice, temps, fields, 10, 42, 1.0);
    
    std::vector<Real> enes(10, 0.5);
    // Need numTypes + 1 vectors (one per type + total)
    std::vector<std::vector<Real>> histMag_x(numTypes + 1, std::vector<Real>(10, 0.0));
    std::vector<std::vector<Real>> histMag_y(numTypes + 1, std::vector<Real>(10, 0.0));
    std::vector<std::vector<Real>> histMag_z(numTypes + 1, std::vector<Real>(10, 1.0));
    
    reporter.partial_report(enes, histMag_x, histMag_y, histMag_z, lattice, 0);
    
    reporter.close();
    
    std::remove("test_partial_output.h5");
    
    std::cout << "PASS\n";
}

void test_reporter_close()
{
    std::cout << "Testing Reporter close... ";
    
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    std::vector<Real> temps = {1.0};
    std::vector<Real> fields = {0.0};
    std::vector<Array> magTypes;
    
    Reporter reporter("test_close_output.h5", magTypes, lattice, temps, fields, 10, 42, 1.0);
    
    reporter.close();
    
    FILE* f = fopen("test_close_output.h5", "r");
    assert(f != nullptr);
    fclose(f);
    
    std::remove("test_close_output.h5");
    
    std::cout << "PASS\n";
}

int main()
{
    std::cout << "Running Reporter tests for vegas...\n\n";
    
    try {
        test_reporter_default_constructor();
        test_reporter_constructor_creates_file();
        test_reporter_file_structure();
        test_reporter_attributes();
        test_reporter_partial_report();
        test_reporter_close();
        
        std::cout << "\nAll Reporter tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}
