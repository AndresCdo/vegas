#ifndef COLOR_BLOCK_H
#define COLOR_BLOCK_H

#include "params.h"
#include <vector>

/**
 * @class ColorBlock
 * @brief Represents a group of atoms with the same color for checkerboard decomposition.
 * 
 * Atoms within a ColorBlock can be updated simultaneously because they are not neighbors
 * (graph coloring guarantees no adjacent atoms share the same color).
 * This enables SIMD vectorization across atoms rather than within neighbor loops.
 */
class ColorBlock {
public:
    ColorBlock() = default;
    
    /**
     * @brief Construct a ColorBlock with a specific color index.
     * @param colorIndex The color index (0, 1, 2, ...)
     */
    explicit ColorBlock(Index colorIndex) : colorIndex_(colorIndex) {}
    
    /**
     * @brief Add an atom to this block.
     * @param atomIndex Global atom index
     */
    void addAtom(Index atomIndex) {
        atomIndices_.push_back(atomIndex);
    }
    
    /**
     * @brief Get the color index of this block.
     */
    Index getColorIndex() const { return colorIndex_; }
    
    /**
     * @brief Get the number of atoms in this block.
     */
    Index getSize() const { return atomIndices_.size(); }
    
    /**
     * @brief Get the global atom index at position i within this block.
     */
    Index getAtomIndex(Index i) const {
        return atomIndices_[i];
    }
    
    /**
     * @brief Get const reference to all atom indices in this block.
     */
    const std::vector<Index>& getAtomIndices() const { return atomIndices_; }
    
    /**
     * @brief Clear all atoms from this block.
     */
    void clear() { atomIndices_.clear(); }
    
    /**
     * @brief Check if block is empty.
     */
    bool empty() const { return atomIndices_.empty(); }
    
private:
    Index colorIndex_ = 0;              ///< Color index (0 = first color)
    std::vector<Index> atomIndices_;    ///< Global atom indices belonging to this block
};

#endif // COLOR_BLOCK_H