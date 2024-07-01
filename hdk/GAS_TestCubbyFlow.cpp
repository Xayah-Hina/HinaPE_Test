#include "GAS_TestCubbyFlow.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_FieldUtils.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <GAS/GAS_Diffuse.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_SparseMatrix.h>

#include <CUDA_CubbyFlow/Core/Solver/Grid/GridSmokeSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Geometry/Box.hpp>
#include <CUDA_CubbyFlow/Core/Geometry/Sphere.hpp>
#include <CUDA_CubbyFlow/Core/Emitter/VolumeGridEmitter3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Advection/CubicSemiLagrangian3.hpp>
#include <CUDA_CubbyFlow/Core/Geometry/RigidBodyCollider.hpp>

#include "Grid/CellCenteredScalarGrid.hpp"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

GAS_TestCubbyFlow::GAS_TestCubbyFlow(const SIM_DataFactory* factory) : BaseClass(factory), Data(nullptr), density_ID(0), temperature_ID(0)
{
}

void GAS_TestCubbyFlow::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
    Data = nullptr;
}

void GAS_TestCubbyFlow::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
    Data = dynamic_cast<const GAS_TestCubbyFlow*>(source)->Data;
}

const SIM_DopDescription* GAS_TestCubbyFlow::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
    ACTIVATE_GAS_DENSITY
    ACTIVATE_GAS_VELOCITY
    ACTIVATE_GAS_SOURCE
    ACTIVATE_GAS_TEMPERATURE
    ACTIVATE_GAS_COLLISION
    ACTIVATE_GAS_PRESSURE
    ACTIVATE_GAS_DIVERGENCE
    PARAMETER_FLOAT(Diffusion, 0.01)
    PARAMETER_FLOAT(ImpulseFactor, 0.1)
    PRMs.emplace_back();

    static SIM_DopDescription DESC(GEN_NODE,
                                   DOP_NAME,
                                   DOP_ENGLISH,
                                   DATANAME,
                                   classname(),
                                   PRMs.data());
    DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
    setGasDescription(DESC);
    return &DESC;
}

bool GAS_TestCubbyFlow::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_ScalarField* D = getScalarField(obj, GAS_NAME_DENSITY);
    SIM_ScalarField* T = getScalarField(obj, GAS_NAME_TEMPERATURE);
    SIM_ScalarField* S = getScalarField(obj, GAS_NAME_SOURCE);
    SIM_VectorField* V = getVectorField(obj, GAS_NAME_VELOCITY);

    if (!D || !T || !V || !S)
    {
        addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
        return false;
    }

    if (!Data)
    {
        using namespace CubbyFlow;
        const Vector3UZ resolution(D->getTotalVoxelRes().x(), D->getTotalVoxelRes().y(), D->getTotalVoxelRes().z());
        const Vector3D spacing(D->getVoxelSize().x(), D->getVoxelSize().y(), D->getVoxelSize().z());
        const Vector3D origin(D->getOrig().x(), D->getOrig().y(), D->getOrig().z());
        Data = std::make_shared<GridSystemData3>();
        Data->Resize(resolution, spacing, origin);

        density_ID = Data->AddAdvectableScalarData(std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);
        temperature_ID = Data->AddAdvectableScalarData(std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);

        std::cout << "resolution: " << resolution.x << " " << resolution.y << " " << resolution.z << '\n';
        std::cout << "spacing: " << spacing.x << " " << spacing.y << " " << spacing.z << '\n';
        std::cout << "origin: " << origin.x << " " << origin.y << " " << origin.z << '\n';
        std::cout << "Data Created!" << '\n';
    }

    return true;
}
