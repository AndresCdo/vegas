#ifndef SPIN_MODEL_TAGS_H
#define SPIN_MODEL_TAGS_H

struct HeisenbergTag {
    static constexpr bool is_ising = false;
    static constexpr bool is_heisenberg = true;
};

struct IsingTag {
    static constexpr bool is_ising = true;
    static constexpr bool is_heisenberg = false;
};

#endif // SPIN_MODEL_TAGS_H
