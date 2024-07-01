#ifndef GAS_TESTCUBBYFLOW_H
#define GAS_TESTCUBBYFLOW_H

#include <CUDA_CubbyFlow/Core/Grid/GridSystemData.hpp>

#include <GAS/GAS_SubSolver.h>

class GAS_TestCubbyFlow : public GAS_SubSolver
{
public:
    inline static const bool GEN_NODE = true;
    inline static const char* DOP_NAME = "TestCubbyFlow";
    inline static const char* DOP_ENGLISH = "Test CubbyFlow";
    inline static const char* DATANAME = "TestCubbyFlow";
    inline static const bool UNIQUE_DATANAME = false;

    std::shared_ptr<CubbyFlow::GridSystemData3> Data;
    size_t density_ID, temperature_ID;

    GETSET_DATA_FUNCS_F("BuoyancyDensity", BuoyancyDensity)
	GETSET_DATA_FUNCS_F("BuoyancyTemperature", BuoyancyTemperature)
	GETSET_DATA_FUNCS_F("Viscosity", Viscosity)

protected:
    explicit GAS_TestCubbyFlow(const SIM_DataFactory* factory);
    void initializeSubclass() final;
    void makeEqualSubclass(const SIM_Data* source) final;
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) final;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_TestCubbyFlow, GAS_SubSolver, "This is a Test CubbyFlow Solver.", getDopDescription());
};

#endif //GAS_TESTCUBBYFLOW_H
