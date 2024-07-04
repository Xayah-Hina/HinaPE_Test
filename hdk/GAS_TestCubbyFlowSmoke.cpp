#include "GAS_TestCubbyFlowSmoke.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_FieldUtils.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <UT/UT_SparseMatrix.h>

#include <CUDA_CubbyFlow/Core/Geometry/Sphere.hpp>
#include <CUDA_CubbyFlow/Core/Emitter/VolumeGridEmitter3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridSmokeSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Advection/CubicSemiLagrangian3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridBackwardEulerDiffusionSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridFractionalSinglePhasePressureSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridBoundaryConditionSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Utils/LevelSetUtils.hpp>
#include <CUDA_CubbyFlow/Core/Array/ArrayUtils.hpp>
#include "CUDA_CubbyFlow/Core/Grid/CellCenteredScalarGrid.hpp"

#include "cubby.h"
#include "Utils/Logging.hpp"

#include <PY/PY_Python.h>
#include <PY/PY_CPythonAPI.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_INT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_INT, 1, &NAME, &Default##NAME);

GAS_TestCubbyFlowSmoke::GAS_TestCubbyFlowSmoke(const SIM_DataFactory* factory) : BaseClass(factory), Data(nullptr), density_ID(0), temperature_ID(0)
{
}

void GAS_TestCubbyFlowSmoke::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
    Data = nullptr;
    density_ID = 0;
    temperature_ID = 0;
}

void GAS_TestCubbyFlowSmoke::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
    Data = dynamic_cast<const GAS_TestCubbyFlowSmoke*>(source)->Data;
    density_ID = dynamic_cast<const GAS_TestCubbyFlowSmoke*>(source)->density_ID;
    temperature_ID = dynamic_cast<const GAS_TestCubbyFlowSmoke*>(source)->temperature_ID;
}

const SIM_DopDescription* GAS_TestCubbyFlowSmoke::getDopDescription()
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
    PARAMETER_FLOAT(BuoyancyDensity, -0.000625)
    PARAMETER_FLOAT(BuoyancyTemperature, 5.0)
    PARAMETER_FLOAT(Viscosity, 0.0)
    PARAMETER_INT(ExtrapolateDepth, 5)
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

bool GAS_TestCubbyFlowSmoke::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
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
    CubbyFlow::Logging::Mute();

    if (!Data)
    {
        const CubbyFlow::Vector3UZ resolution(D->getTotalVoxelRes().x(), D->getTotalVoxelRes().y(), D->getTotalVoxelRes().z());
        const CubbyFlow::Vector3D spacing(D->getVoxelSize().x(), D->getVoxelSize().y(), D->getVoxelSize().z());
        const CubbyFlow::Vector3D origin(D->getOrig().x(), D->getOrig().y(), D->getOrig().z());
        Data = std::make_shared<CubbyFlow::GridSystemData3>();
        Data->Resize(resolution, spacing, origin);

        density_ID = Data->AddAdvectableScalarData(std::make_shared<CubbyFlow::CellCenteredScalarGrid3::Builder>(), 0.0);
        temperature_ID = Data->AddAdvectableScalarData(std::make_shared<CubbyFlow::CellCenteredScalarGrid3::Builder>(), 0.0);
    }

    static std::shared_ptr<CubbyFlow::CubicSemiLagrangian3> AdvectionSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridBackwardEulerDiffusionSolver3> DiffusionSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridFractionalSinglePhasePressureSolver3> PressureSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridBoundaryConditionSolver3> BoundaryConditionSolver = nullptr;
    if (!AdvectionSolver) AdvectionSolver = std::make_shared<CubbyFlow::CubicSemiLagrangian3>();
    if (!DiffusionSolver) DiffusionSolver = std::make_shared<CubbyFlow::GridBackwardEulerDiffusionSolver3>();
    if (!PressureSolver) PressureSolver = std::make_shared<CubbyFlow::GridFractionalSinglePhasePressureSolver3>();
    if (!BoundaryConditionSolver) BoundaryConditionSolver = PressureSolver->SuggestedBoundaryConditionSolver();

    CubbyFlow::ScalarGrid3Ptr D_Cubby = Data->AdvectableScalarDataAt(density_ID);
    CubbyFlow::ScalarGrid3Ptr T_Cubby = Data->AdvectableScalarDataAt(temperature_ID);
    CubbyFlow::FaceCenteredGrid3Ptr V_Cubby = Data->Velocity();
    CubbyFlow::ScalarField3Ptr FluidSDF_Cubby = std::make_shared<CubbyFlow::ConstantScalarField3>(-std::numeric_limits<double>::max());
    const auto depth = getExtrapolateDepth();


    // Emit Source
    CubbyFlow::EmitFromHoudiniSource(D_Cubby, S);
    CubbyFlow::EmitFromHoudiniSource(T_Cubby, S);


    // Buoyancy Force
    {
        double tAmb = 0.0;
        T_Cubby->ForEachCellIndex([&](size_t i, size_t j, size_t k) { tAmb += (*T_Cubby)(i, j, k); });
        tAmb /= static_cast<double>(T_Cubby->Resolution().x * T_Cubby->Resolution().y * T_Cubby->Resolution().z);

        CubbyFlow::ArrayView3<double> u = V_Cubby->UView();
        CubbyFlow::ArrayView3<double> v = V_Cubby->VView();
        CubbyFlow::ArrayView3<double> w = V_Cubby->WView();
        auto uPos = V_Cubby->UPosition();
        auto vPos = V_Cubby->VPosition();
        auto wPos = V_Cubby->WPosition();

        V_Cubby->ParallelForEachVIndex([&](const CubbyFlow::Vector3UZ& idx)
        {
            const CubbyFlow::Vector3D pt = vPos(idx);
            const double fBuoy = getBuoyancyDensity() * D_Cubby->Sample(pt) + getBuoyancyTemperature() * (T_Cubby->Sample(pt) - tAmb);
            v(idx) += timestep * fBuoy;
        });

        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);
    }


    // Viscosity
    {
        const std::shared_ptr<CubbyFlow::FaceCenteredGrid3> vel0 = std::dynamic_pointer_cast<CubbyFlow::FaceCenteredGrid3>(V_Cubby->Clone());
        DiffusionSolver->Solve(*vel0, getViscosity(), timestep, V_Cubby.get(), *BoundaryConditionSolver->GetColliderSDF(), *FluidSDF_Cubby);
        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);
    }


    // Pressure
    {
        const std::shared_ptr<CubbyFlow::FaceCenteredGrid3> vel0 = std::dynamic_pointer_cast<CubbyFlow::FaceCenteredGrid3>(V_Cubby->Clone());
        PressureSolver->Solve(*vel0, timestep, V_Cubby.get(), *BoundaryConditionSolver->GetColliderSDF(), *BoundaryConditionSolver->GetColliderVelocityField(), *FluidSDF_Cubby, false);
        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);
    }


    // Advection
    {
        size_t n = Data->NumberOfAdvectableScalarData();
        for (size_t i = 0; i < n; ++i)
        {
            CubbyFlow::ScalarGrid3Ptr grid = Data->AdvectableScalarDataAt(i);
            std::shared_ptr<CubbyFlow::ScalarGrid3> grid0 = grid->Clone();

            AdvectionSolver->Advect(*grid0, *V_Cubby, timestep, grid.get(), *BoundaryConditionSolver->GetColliderSDF());

            CubbyFlow::Array3<char> marker(grid->DataSize());
            CubbyFlow::GridDataPositionFunc<3> pos = grid->DataPosition();
            ParallelForEachIndex(marker.Size(), [&](size_t i, size_t j, size_t k)
            {
                if (CubbyFlow::IsInsideSDF(BoundaryConditionSolver->GetColliderSDF()->Sample(pos(i, j, k))))
                {
                    marker(i, j, k) = 0;
                } else
                {
                    marker(i, j, k) = 1;
                }
            });
            CubbyFlow::ExtrapolateToRegion(grid->DataView(), marker, depth, grid->DataView());
        }

        n = Data->NumberOfAdvectableVectorData();
        const std::shared_ptr<CubbyFlow::FaceCenteredGrid3> vel0 = std::dynamic_pointer_cast<CubbyFlow::FaceCenteredGrid3>(V_Cubby->Clone());
        AdvectionSolver->Advect(*vel0, *vel0, timestep, V_Cubby.get(), *BoundaryConditionSolver->GetColliderSDF());
        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);
    }

    CubbyFlow::WriteHoudiniField(D, D_Cubby);
    CubbyFlow::WriteHoudiniField(T, T_Cubby);
    CubbyFlow::WriteHoudiniField(V, V_Cubby);

    UT_WorkBuffer expr;
    PY_Result result;

    expr.sprintf(R"(
from phi.torch.flow import *

velocity = StaggeredGrid((0, 0, 0), 0, x=50, y=50, z=50, bounds=Box(x=1, y=1, z=1))  # or CenteredGrid(...)
smoke = CenteredGrid(0, ZERO_GRADIENT, x=50, y=50, z=50, bounds=Box(x=1, y=1, z=1))
INFLOW = 0.5 * resample(Sphere(x=0.5, y=0.5, z=0.5, radius=0.1), to=smoke, soft=True)
pressure = None


@jit_compile  # Only for PyTorch, TensorFlow and Jax
def step(v, s, p, dt=1.):
    s = advect.mac_cormack(s, v, dt) + INFLOW
    buoyancy = resample(s * (0, 0, 0.1), to=v)
    v = advect.semi_lagrangian(v, v, dt) + buoyancy * dt
    v, p = fluid.make_incompressible(v, (), Solve('auto', 1e-5, x0=p))
    return v, s, p


velocity, smoke, pressure = step(velocity, smoke, pressure)
d_n = smoke.data.native('x,y,z').cpu().numpy().flatten()
d_n_no_zero = d_n[d_n != 0]
)");
    PYrunPythonStatementsAndExpectNoErrors(expr.buffer());
    expr.sprintf("d_n_no_zero.tolist()\n");
    result = PYrunPythonExpressionAndExpectNoErrors(expr.buffer(), PY_Result::DOUBLE_ARRAY);
    if (result.myResultType != PY_Result::DOUBLE_ARRAY)
    {
        printf("Error: %s\n", result.myErrValue.buffer());
        return false;
    }

    printf("Python Result:");
    for (double i : result.myDoubleArray) printf(" %f", i);
    printf("\n");

    return true;
}
