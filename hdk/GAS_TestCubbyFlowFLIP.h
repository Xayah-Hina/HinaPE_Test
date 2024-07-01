#ifndef GAS_TESTCUBBYFLOWFLIP_H
#define GAS_TESTCUBBYFLOWFLIP_H

#include <CUDA_CubbyFlow/Core/Particle/ParticleSystemData.hpp>
#include <CUDA_CubbyFlow/Core/Grid/GridSystemData.hpp>

#include <GAS/GAS_SubSolver.h>

class GAS_TestCubbyFlowFLIP : public GAS_SubSolver
{
public:
    inline static const bool GEN_NODE = true;
    inline static const char* DOP_NAME = "TestCubbyFlowFLIP";
    inline static const char* DOP_ENGLISH = "Test CubbyFlow FLIP";
    inline static const char* DATANAME = "TestCubbyFlowFLIP";
    inline static const bool UNIQUE_DATANAME = false;

    std::shared_ptr<CubbyFlow::ParticleSystemData3> Particles;
    std::shared_ptr<CubbyFlow::GridSystemData3> Data;
    size_t sdf_ID;

    GETSET_DATA_FUNCS_F("Gravity", Gravity)
    GETSET_DATA_FUNCS_F("Viscosity", Viscosity)
	GETSET_DATA_FUNCS_I("ExtrapolateDepth", ExtrapolateDepth)

protected:
    explicit GAS_TestCubbyFlowFLIP(const SIM_DataFactory* factory);
    void initializeSubclass() final;
    void makeEqualSubclass(const SIM_Data* source) final;
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) final;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_TestCubbyFlowFLIP, GAS_SubSolver, "This is a Test CubbyFlow FLIP Solver.", getDopDescription());
};


#endif //GAS_TESTCUBBYFLOWFLIP_H
