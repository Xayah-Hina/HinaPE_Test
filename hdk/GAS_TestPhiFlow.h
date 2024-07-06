#ifndef GAS_TESTPHIFLOW_H
#define GAS_TESTPHIFLOW_H

#include <GAS/GAS_SubSolver.h>

class GAS_TestPhiFlow : public GAS_SubSolver
{
public:
    inline static const bool GEN_NODE = true;
    inline static const char* DOP_NAME = "TestPhiFlow";
    inline static const char* DOP_ENGLISH = "Test PhiFlow";
    inline static const char* DATANAME = "TestPhiFlow";
    inline static const bool UNIQUE_DATANAME = false;

protected:
    explicit GAS_TestPhiFlow(const SIM_DataFactory* factory);
    void initializeSubclass() final;
    void makeEqualSubclass(const SIM_Data* source) final;
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) final;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_TestPhiFlow, GAS_SubSolver, "This is a Test PhiFlow Smoke Solver.", getDopDescription());
};

#endif //GAS_TESTPHIFLOW_H
