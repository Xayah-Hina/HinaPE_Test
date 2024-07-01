#ifndef HINAPE_TEST_CUBBY_H
#define HINAPE_TEST_CUBBY_H

#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <GU/GU_Detail.h>

#include <CUDA_CubbyFlow/Core/Grid/ScalarGrid.hpp>
#include <CUDA_CubbyFlow/Core/Grid/FaceCenteredGrid.hpp>
#include <CUDA_CubbyFlow/Core/Particle/ParticleSystemData.hpp>

namespace CubbyFlow
{
void EmitFromHoudiniSource(std::shared_ptr<ScalarGrid3> &TARGET, const SIM_ScalarField *SOURCE);
void EmitFromHoudiniSource(std::shared_ptr<FaceCenteredGrid3> &TARGET, const SIM_VectorField *SOURCE);
void WriteHoudiniField(SIM_ScalarField *TARGET, const std::shared_ptr<ScalarGrid3> &SOURCE);
void WriteHoudiniField(SIM_VectorField *TARGET, const std::shared_ptr<FaceCenteredGrid3> &SOURCE);

void LoadHoudiniParticles(std::shared_ptr<CubbyFlow::ParticleSystemData3>& Particles, const GU_Detail& gdp);
void WriteHoudiniParticles(GU_Detail& gdp, const std::shared_ptr<CubbyFlow::ParticleSystemData3>& Particles);
}

#endif //HINAPE_TEST_CUBBY_H
