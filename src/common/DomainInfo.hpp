#pragma once

#include "common/Vec3.hpp"

// Simulation dimension, independent of field storage or rendering format.
enum class Dimension {
    D2 = 2,
    D3 = 3
};

struct DomainBounds {
    Vec3 min;
    Vec3 max;
};

