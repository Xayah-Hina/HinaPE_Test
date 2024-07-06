#include "GAS_TestPhiFlow.h"

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


#include <PY/PY_Python.h>
#include <PY/PY_CPythonAPI.h>

#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);

static int To1DIdx(const UT_Vector3I& idx, const UT_Vector3I& res)
{
    return idx.x() + res.x() * (idx.y() + res.y() * idx.z());
}

GAS_TestPhiFlow::GAS_TestPhiFlow(const SIM_DataFactory* factory): BaseClass(factory)
{
}

void GAS_TestPhiFlow::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();

    UT_WorkBuffer expr;
    expr.sprintf(R"(
from phi.torch.flow import *
res = (100, 100, 100)
velocity = StaggeredGrid((0, 0, 0), 0, x=res[0], y=res[1], z=res[2], bounds=Box(x=1, y=1, z=1))
smoke = CenteredGrid(0, ZERO_GRADIENT, x=res[0], y=res[1], z=res[2], bounds=Box(x=1, y=1, z=1))
INFLOW = resample(Sphere(x=0.5, y=0.5, z=0.5, radius=0.1), to=smoke, soft=True)
pressure = None

@jit_compile  # Only for PyTorch, TensorFlow and Jax
def step(v, s, p, dt=1.):
    s = advect.mac_cormack(s, v, dt) + INFLOW
    buoyancy = resample(s * (0, 0, 0.1), to=v)
    v = advect.semi_lagrangian(v, v, dt) + buoyancy * dt
    v, p = fluid.make_incompressible(v, (), Solve('auto', 1e-5, x0=p))
    return v, s, p
)");
    PYrunPythonStatementsAndExpectNoErrors(expr.buffer());
}

void GAS_TestPhiFlow::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
}

const SIM_DopDescription* GAS_TestPhiFlow::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_DENSITY
    ACTIVATE_GAS_VELOCITY
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

void WriteFieldPartial(SIM_RawField* TARGET, const std::vector<double>& SOURCE, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setArray(TARGET->fieldNC());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        UT_Vector3I cell(vit.x(), vit.y(), vit.z());
        int idx = To1DIdx(cell, TARGET->getVoxelRes());
        vit.setValue(static_cast<float>(SOURCE[idx]));
    }
}

THREADED_METHOD2(, true, WriteField, SIM_RawField *, TARGET, const std::vector<double>&, SOURCE);

void WriteHoudiniField(SIM_ScalarField* TARGET, const std::vector<double>& SOURCE)
{
    const auto f = TARGET->getField()->field();
    auto voxels = f->getXRes() * f->getYRes() * f->getZRes();
    if (voxels != SOURCE.size())
    {
        printf("Error: Field size mismatch\n");
        return;
    }
    WriteField(TARGET->getField(), SOURCE);
}


bool GAS_TestPhiFlow::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_ScalarField* D = getScalarField(obj, GAS_NAME_DENSITY);
    SIM_VectorField* V = getVectorField(obj, GAS_NAME_VELOCITY);
    if (!D || !V)
    {
        addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
        return false;
    }

    UT_WorkBuffer expr;
    PY_Result result;
    expr.sprintf(R"(
velocity, smoke, pressure = step(velocity, smoke, pressure)
smoke_output = smoke.data.native('x,z,y').cpu().numpy().flatten()
)");
    PYrunPythonStatementsAndExpectNoErrors(expr.buffer());
    expr.sprintf("smoke_output.tolist()\n");
    result = PYrunPythonExpressionAndExpectNoErrors(expr.buffer(), PY_Result::DOUBLE_ARRAY);
    if (result.myResultType != PY_Result::DOUBLE_ARRAY)
    {
        printf("Error: %s\n", result.myErrValue.buffer());
        return false;
    }

    WriteHoudiniField(D, result.myDoubleArray);

    return true;
}
