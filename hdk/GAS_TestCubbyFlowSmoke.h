#ifndef GAS_TESTCUBBYFLOW_SMOKE_H
#define GAS_TESTCUBBYFLOW_SMOKE_H

#include <CUDA_CubbyFlow/Core/Grid/GridSystemData.hpp>

#include <GAS/GAS_SubSolver.h>

class GAS_TestCubbyFlowSmoke : public GAS_SubSolver
{
public:
    inline static const bool GEN_NODE = true;
    inline static const char* DOP_NAME = "TestCubbyFlowSmoke";
    inline static const char* DOP_ENGLISH = "Test CubbyFlow Smoke";
    inline static const char* DATANAME = "TestCubbyFlowSmoke";
    inline static const bool UNIQUE_DATANAME = false;

    std::shared_ptr<CubbyFlow::GridSystemData3> Data;
    size_t density_ID, temperature_ID;

    GETSET_DATA_FUNCS_F("BuoyancyDensity", BuoyancyDensity)
	GETSET_DATA_FUNCS_F("BuoyancyTemperature", BuoyancyTemperature)
	GETSET_DATA_FUNCS_F("Viscosity", Viscosity)
	GETSET_DATA_FUNCS_I("ExtrapolateDepth", ExtrapolateDepth)

protected:
    explicit GAS_TestCubbyFlowSmoke(const SIM_DataFactory* factory);
    void initializeSubclass() final;
    void makeEqualSubclass(const SIM_Data* source) final;
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) final;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_TestCubbyFlowSmoke, GAS_SubSolver, "This is a Test CubbyFlow Smoke Solver.", getDopDescription());
};

#endif //GAS_TESTCUBBYFLOW_SMOKE_H
