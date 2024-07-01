#ifndef HINAPE_TEST_CUBBY_H
#define HINAPE_TEST_CUBBY_H

#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>

#include <CUDA_CubbyFlow/Core/Grid/ScalarGrid.hpp>
#include <CUDA_CubbyFlow/Core/Grid/FaceCenteredGrid.hpp>

#include "Grid/FaceCenteredGrid.hpp"

namespace CubbyFlow
{
void EmitFromHoudiniSource(std::shared_ptr<ScalarGrid3> &TARGET, const SIM_ScalarField *SOURCE);
void EmitFromHoudiniSource(std::shared_ptr<FaceCenteredGrid3> &TARGET, const SIM_VectorField *SOURCE);
void WriteHoudiniField(SIM_ScalarField *TARGET, const std::shared_ptr<ScalarGrid3> &SOURCE);
void WriteHoudiniField(SIM_VectorField *TARGET, const std::shared_ptr<FaceCenteredGrid3> &SOURCE);
}

#endif //HINAPE_TEST_CUBBY_H
