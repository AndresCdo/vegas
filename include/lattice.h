#ifndef LATTICE
#define LATTICE

#include <vector>
#include <fstream>
#include <map>

#include "atom.h"

class Lattice
{
public:
    Lattice(std::string fileName);
    ~Lattice();
    
    // Move constructor and assignment
    Lattice(Lattice&& other) noexcept;
    Lattice& operator=(Lattice&& other) noexcept;
    
    // Delete copy constructor and assignment
    Lattice(const Lattice&) = delete;
    Lattice& operator=(const Lattice&) = delete;

    std::vector<Atom>& getAtoms();

    const std::map<std::string, Index>& getMapTypeIndexes() const;
    const std::map<Index, std::string>& getMapIndexTypes() const;
    const std::vector<Index>& getSizesByIndex() const;

private:
    std::vector<Atom> atoms_;
    std::map<std::string, Index> mapTypeIndexes_;
    std::map<Index, std::string> mapIndexTypes_;
    std::vector<Index> sizesByIndex_;

};

#endif