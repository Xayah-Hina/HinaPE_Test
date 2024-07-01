#include "GAS_TestCubbyFlowFLIP.h"

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

#include <CUDA_CubbyFlow/Core/Grid/CellCenteredScalarGrid.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Advection/CubicSemiLagrangian3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridBackwardEulerDiffusionSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridFractionalSinglePhasePressureSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Solver/Grid/GridBoundaryConditionSolver3.hpp>
#include <CUDA_CubbyFlow/Core/Utils/LevelSetUtils.hpp>
#include <CUDA_CubbyFlow/Core/Array/ArrayUtils.hpp>

#include "cubby.h"
#include "Utils/Logging.hpp"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);
#define ACTIVATE_GAS_SURFACE static PRM_Name SurfaceName(GAS_NAME_SURFACE, "Surface"); static PRM_Default SurfaceNameDefault(0, GAS_NAME_SURFACE); PRMs.emplace_back(PRM_STRING, 1, &SurfaceName, &SurfaceNameDefault);

#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_INT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_INT, 1, &NAME, &Default##NAME);

GAS_TestCubbyFlowFLIP::GAS_TestCubbyFlowFLIP(const SIM_DataFactory* factory) : BaseClass(factory), Particles(nullptr), Data(nullptr), sdf_ID(0)
{
}

void GAS_TestCubbyFlowFLIP::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
    Particles = nullptr;
    Data = nullptr;
    sdf_ID = 0;
}

void GAS_TestCubbyFlowFLIP::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
    Particles = dynamic_cast<const GAS_TestCubbyFlowFLIP*>(source)->Particles;
    Data = dynamic_cast<const GAS_TestCubbyFlowFLIP*>(source)->Data;
    sdf_ID = dynamic_cast<const GAS_TestCubbyFlowFLIP*>(source)->sdf_ID;
}

const SIM_DopDescription* GAS_TestCubbyFlowFLIP::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
    ACTIVATE_GAS_VELOCITY
    ACTIVATE_GAS_COLLISION
    ACTIVATE_GAS_PRESSURE
    ACTIVATE_GAS_DIVERGENCE
    ACTIVATE_GAS_SURFACE
    PARAMETER_FLOAT(Gravity, -9.8)
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

bool GAS_TestCubbyFlowFLIP::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_VectorField* V = getVectorField(obj, GAS_NAME_VELOCITY);
    SIM_ScalarField* SF = getScalarField(obj, GAS_NAME_SURFACE);

    if (!G || !V || !SF)
    {
        addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
        return false;
    }
    CubbyFlow::Logging::Mute();

    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    if (gdp.getNumPoints() == 0) return true;

    if (!Data)
    {
        const CubbyFlow::Vector3UZ resolution(V->getTotalVoxelRes().x(), V->getTotalVoxelRes().y(), V->getTotalVoxelRes().z());
        const CubbyFlow::Vector3D spacing(V->getVoxelSize().x(), V->getVoxelSize().y(), V->getVoxelSize().z());
        const CubbyFlow::Vector3D origin(V->getOrig().x(), V->getOrig().y(), V->getOrig().z());
        Particles = std::make_shared<CubbyFlow::ParticleSystemData3>();
        Data = std::make_shared<CubbyFlow::GridSystemData3>();
        Data->Resize(resolution, spacing, origin);

        sdf_ID = Data->AddAdvectableScalarData(std::make_shared<CubbyFlow::CellCenteredScalarGrid3::Builder>(), 0.0);

        // Load Particles Once
        CubbyFlow::LoadHoudiniParticles(Particles, gdp);
    }

    static std::shared_ptr<CubbyFlow::CubicSemiLagrangian3> AdvectionSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridBackwardEulerDiffusionSolver3> DiffusionSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridFractionalSinglePhasePressureSolver3> PressureSolver = nullptr;
    static std::shared_ptr<CubbyFlow::GridBoundaryConditionSolver3> BoundaryConditionSolver = nullptr;
    if (!AdvectionSolver) AdvectionSolver = std::make_shared<CubbyFlow::CubicSemiLagrangian3>();
    if (!DiffusionSolver) DiffusionSolver = std::make_shared<CubbyFlow::GridBackwardEulerDiffusionSolver3>();
    if (!PressureSolver) PressureSolver = std::make_shared<CubbyFlow::GridFractionalSinglePhasePressureSolver3>();
    if (!BoundaryConditionSolver) BoundaryConditionSolver = PressureSolver->SuggestedBoundaryConditionSolver();

    CubbyFlow::ScalarGrid3Ptr FluidSDF_Cubby = Data->AdvectableScalarDataAt(sdf_ID);
    CubbyFlow::FaceCenteredGrid3Ptr V_Cubby = Data->Velocity();
    CubbyFlow::ArrayView1<CubbyFlow::Vector3<double>> positions = Particles->Positions();
    CubbyFlow::ArrayView1<CubbyFlow::Vector3<double>> velocities = Particles->Velocities();
    CubbyFlow::Array3<char> m_uMarkers;
    CubbyFlow::Array3<char> m_vMarkers;
    CubbyFlow::Array3<char> m_wMarkers;
    CubbyFlow::Array3<double> m_uDelta;
    CubbyFlow::Array3<double> m_vDelta;
    CubbyFlow::Array3<double> m_wDelta;
    const size_t N = Particles->NumberOfParticles();
    const auto depth = getExtrapolateDepth();

    BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth); // Update All First

    // P2G
    {
        V_Cubby->Fill(CubbyFlow::Vector3D{});
        CubbyFlow::ArrayView3<double> u = V_Cubby->UView();
        CubbyFlow::ArrayView3<double> v = V_Cubby->VView();
        CubbyFlow::ArrayView3<double> w = V_Cubby->WView();
        CubbyFlow::Array3<double> uWeight{u.Size()};
        CubbyFlow::Array3<double> vWeight{v.Size()};
        CubbyFlow::Array3<double> wWeight{w.Size()};
        m_uMarkers.Resize(u.Size());
        m_vMarkers.Resize(v.Size());
        m_wMarkers.Resize(w.Size());
        m_uMarkers.Fill(0);
        m_vMarkers.Fill(0);
        m_wMarkers.Fill(0);

        CubbyFlow::LinearArraySampler3<double> uSampler{V_Cubby->UView(), V_Cubby->GridSpacing(), V_Cubby->UOrigin()};
        CubbyFlow::LinearArraySampler3<double> vSampler{V_Cubby->VView(), V_Cubby->GridSpacing(), V_Cubby->VOrigin()};
        CubbyFlow::LinearArraySampler3<double> wSampler{V_Cubby->WView(), V_Cubby->GridSpacing(), V_Cubby->WOrigin()};

        for (size_t i = 0; i < N; ++i)
        {
            std::array<CubbyFlow::Vector3UZ, 8> indices{};
            std::array<double, 8> weights{};

            uSampler.GetCoordinatesAndWeights(positions[i], indices, weights);
            for (int j = 0; j < 8; ++j)
            {
                u(indices[j]) += velocities[i].x * weights[j];
                uWeight(indices[j]) += weights[j];
                m_uMarkers(indices[j]) = 1;
            }

            vSampler.GetCoordinatesAndWeights(positions[i], indices, weights);
            for (int j = 0; j < 8; ++j)
            {
                v(indices[j]) += velocities[i].y * weights[j];
                vWeight(indices[j]) += weights[j];
                m_vMarkers(indices[j]) = 1;
            }

            wSampler.GetCoordinatesAndWeights(positions[i], indices, weights);
            for (int j = 0; j < 8; ++j)
            {
                w(indices[j]) += velocities[i].z * weights[j];
                wWeight(indices[j]) += weights[j];
                m_wMarkers(indices[j]) = 1;
            }
        }

        ParallelForEachIndex(uWeight.Size(), [&](size_t i, size_t j, size_t k)
        {
            if (uWeight(i, j, k) > 0.0)
            {
                u(i, j, k) /= uWeight(i, j, k);
            }
        });
        ParallelForEachIndex(vWeight.Size(), [&](size_t i, size_t j, size_t k)
        {
            if (vWeight(i, j, k) > 0.0)
            {
                v(i, j, k) /= vWeight(i, j, k);
            }
        });
        ParallelForEachIndex(wWeight.Size(), [&](size_t i, size_t j, size_t k)
        {
            if (wWeight(i, j, k) > 0.0)
            {
                w(i, j, k) /= wWeight(i, j, k);
            }
        });

        // Store snapshot
        m_uDelta.Resize(u.Size());
        m_vDelta.Resize(v.Size());
        m_wDelta.Resize(w.Size());
        V_Cubby->ParallelForEachUIndex([&](const CubbyFlow::Vector3UZ& idx) { m_uDelta(idx) = u(idx); });
        V_Cubby->ParallelForEachVIndex([&](const CubbyFlow::Vector3UZ& idx) { m_vDelta(idx) = v(idx); });
        V_Cubby->ParallelForEachWIndex([&](const CubbyFlow::Vector3UZ& idx) { m_wDelta(idx) = w(idx); });
    }


    // Build SDF
    {
        CubbyFlow::GridDataPositionFunc<3> sdfPos = FluidSDF_Cubby->DataPosition();
        const double maxH = std::max({FluidSDF_Cubby->GridSpacing().x, FluidSDF_Cubby->GridSpacing().y, FluidSDF_Cubby->GridSpacing().z});
        double radius = 1.2 * maxH / std::sqrt(2.0);
        double sdfBandRadius = 2.0 * radius;
        Particles->BuildNeighborSearcher(2 * radius);
        CubbyFlow::PointNeighborSearcher3Ptr searcher = Particles->NeighborSearcher();
        FluidSDF_Cubby->ParallelForEachDataPointIndex([&](size_t i, size_t j, size_t k)
        {
            CubbyFlow::Vector3D pt = sdfPos(i, j, k);
            double minDist = sdfBandRadius;

            searcher->ForEachNearbyPoint(
                pt, sdfBandRadius, [&](size_t, const CubbyFlow::Vector3D& x)
                {
                    minDist = std::min(minDist, pt.DistanceTo(x));
                });
            (*FluidSDF_Cubby)(i, j, k) = minDist - radius;
        });

        CubbyFlow::Array3<char> marker(FluidSDF_Cubby->DataSize());
        CubbyFlow::GridDataPositionFunc<3> pos = FluidSDF_Cubby->DataPosition();
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
        CubbyFlow::ExtrapolateToRegion(FluidSDF_Cubby->DataView(), marker, depth, FluidSDF_Cubby->DataView());
    }


    // Extrapolate to Air
    {
        const CubbyFlow::ArrayView3<double> u = V_Cubby->UView();
        const CubbyFlow::ArrayView3<double> v = V_Cubby->VView();
        const CubbyFlow::ArrayView3<double> w = V_Cubby->WView();

        ExtrapolateToRegion(V_Cubby->UView(), m_uMarkers, depth, u);
        ExtrapolateToRegion(V_Cubby->VView(), m_vMarkers, depth, v);
        ExtrapolateToRegion(V_Cubby->WView(), m_wMarkers, depth, w);

        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);
    }


    // Gravity Force
    {
        V_Cubby->ForEachVIndex([&](const CubbyFlow::Vector3UZ& idx)
        {
            V_Cubby->VView()(idx) += timestep * getGravity();
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
        const CubbyFlow::ArrayView3<double> u = V_Cubby->UView();
        const CubbyFlow::ArrayView3<double> v = V_Cubby->VView();
        const CubbyFlow::ArrayView3<double> w = V_Cubby->WView();

        ExtrapolateToRegion(V_Cubby->UView(), m_uMarkers, depth, u);
        ExtrapolateToRegion(V_Cubby->VView(), m_vMarkers, depth, v);
        ExtrapolateToRegion(V_Cubby->WView(), m_wMarkers, depth, w);

        BoundaryConditionSolver->ConstrainVelocity(V_Cubby.get(), depth);

        // G2P
        V_Cubby->ParallelForEachUIndex([&](const CubbyFlow::Vector3UZ& idx) { m_uDelta(idx) = V_Cubby->U(idx) - m_uDelta(idx); });
        V_Cubby->ParallelForEachVIndex([&](const CubbyFlow::Vector3UZ& idx) { m_vDelta(idx) = V_Cubby->V(idx) - m_vDelta(idx); });
        V_Cubby->ParallelForEachWIndex([&](const CubbyFlow::Vector3UZ& idx) { m_wDelta(idx) = V_Cubby->W(idx) - m_wDelta(idx); });
        CubbyFlow::LinearArraySampler3<double> uSampler{m_uDelta.View(), V_Cubby->GridSpacing().CastTo<double>(), V_Cubby->UOrigin().CastTo<double>()};
        CubbyFlow::LinearArraySampler3<double> vSampler{m_vDelta.View(), V_Cubby->GridSpacing().CastTo<double>(), V_Cubby->VOrigin().CastTo<double>()};
        CubbyFlow::LinearArraySampler3<double> wSampler{m_wDelta.View(), V_Cubby->GridSpacing().CastTo<double>(), V_Cubby->WOrigin().CastTo<double>()};

        auto sampler = [uSampler, vSampler, wSampler](const CubbyFlow::Vector3D& x)
        {
            const CubbyFlow::Vector3<double> xf = x.CastTo<double>();
            const double u = uSampler(xf);
            const double v = vSampler(xf);
            const double w = wSampler(xf);
            return CubbyFlow::Vector3D{u, v, w};
        };

        CubbyFlow::ParallelFor(CubbyFlow::ZERO_SIZE, N, [&](size_t i)
        {
            CubbyFlow::Vector3D flipVel = velocities[i] + sampler(positions[i]);
            velocities[i] = flipVel;
        });


        CubbyFlow::BoundingBox3D boundingBox = V_Cubby->GetBoundingBox();
        CubbyFlow::ParallelFor(CubbyFlow::ZERO_SIZE, N, [&](size_t i)
        {
            CubbyFlow::Vector3D pt0 = positions[i];
            CubbyFlow::Vector3D pt1 = pt0;
            CubbyFlow::Vector3D vel = velocities[i];

            // Adaptive time-stepping
            const unsigned int numSubSteps = static_cast<unsigned int>(std::max(5., 1.0));
            const double dt = timestep / numSubSteps;
            for (unsigned int t = 0; t < numSubSteps; ++t)
            {
                CubbyFlow::Vector3D vel0 = V_Cubby->Sample(pt0);

                // Mid-point rule
                CubbyFlow::Vector3D midPt = pt0 + 0.5 * dt * vel0;
                CubbyFlow::Vector3D midVel = V_Cubby->Sample(midPt);
                pt1 = pt0 + dt * midVel;

                pt0 = pt1;
            }

            if (pt1.x <= boundingBox.lowerCorner.x)
            {
                pt1.x = boundingBox.lowerCorner.x;
                vel.x = 0.0;
            }
            if (pt1.x >= boundingBox.upperCorner.x)
            {
                pt1.x = boundingBox.upperCorner.x;
                vel.x = 0.0;
            }
            if (pt1.y <= boundingBox.lowerCorner.y)
            {
                pt1.y = boundingBox.lowerCorner.y;
                vel.y = 0.0;
            }
            if (pt1.y >= boundingBox.upperCorner.y)
            {
                pt1.y = boundingBox.upperCorner.y;
                vel.y = 0.0;
            }
            if (pt1.z <= boundingBox.lowerCorner.z)
            {
                pt1.z = boundingBox.lowerCorner.z;
                vel.z = 0.0;
            }
            if (pt1.z >= boundingBox.upperCorner.z)
            {
                pt1.z = boundingBox.upperCorner.z;
                vel.z = 0.0;
            }

            positions[i] = pt1;
            velocities[i] = vel;
        });
    }

    CubbyFlow::WriteHoudiniParticles(gdp, Particles);
    CubbyFlow::WriteHoudiniField(SF, FluidSDF_Cubby);
    CubbyFlow::WriteHoudiniField(V, V_Cubby);

    return true;
}
