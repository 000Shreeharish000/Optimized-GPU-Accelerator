// Crossover from an interior / first-order primal-dual point to a simplex basis.
//
// Basis identification: every column and row gets an "interior-ness" score
// (distance to the nearest bound relative to the magnitude of its dual); the m
// highest-scoring entities become basic, the others nonbasic at their nearest
// bound. The simplex (lp/simplex.h) then repairs singularities (logical
// substitution) and removes the remaining primal/dual infeasibilities from this
// warm start. For a near-optimal input point this needs few pivots.
#pragma once

#include <vector>

#include "core/model.h"
#include "lp/simplex.h"
#include "util/timer.h"

namespace pramana {

void guessBasis(const Model& model, const std::vector<double>& x, const std::vector<double>& rowDual,
                std::vector<BasisStatus>& colStatus, std::vector<BasisStatus>& rowStatus);

// Guess a basis from (x, y) and finish with the simplex. Returns the simplex result.
LpResult crossover(const Model& model, const std::vector<double>& x, const std::vector<double>& rowDual,
                   const SimplexOptions& opts, const Deadline* deadline);

}  // namespace pramana
