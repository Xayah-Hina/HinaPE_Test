#ifndef GAS_READREFIMAGE_H
#define GAS_READREFIMAGE_H


#include <GAS/GAS_SubSolver.h>

class GAS_ReadRefImage : public GAS_SubSolver
{
public:
    inline static const bool GEN_NODE = true;
    inline static const char* DOP_NAME = "ReadRefImage";
    inline static const char* DOP_ENGLISH = "Read Ref Image";
    inline static const char* DATANAME = "ReadRefImage";
    inline static const bool UNIQUE_DATANAME = false;

    GET_DATA_FUNC_S("RefImage", RefImage)

protected:
    explicit GAS_ReadRefImage(const SIM_DataFactory* factory);
    void initializeSubclass() final;
    void makeEqualSubclass(const SIM_Data* source) final;
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) final;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_ReadRefImage, GAS_SubSolver, "This is a Read Reference Image Solver.", getDopDescription());
};

#endif //GAS_READREFIMAGE_H
