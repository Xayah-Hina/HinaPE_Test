#include "GAS_ReadRefImage.h"

#include <SIM/SIM_Engine.h>
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

#define GAS_NAME_REFERENCE		"reference"
#define ACTIVATE_GAS_REFERENCE static PRM_Name ReferenceName(GAS_NAME_REFERENCE, "Reference"); static PRM_Default ReferenceNameDefault(0, GAS_NAME_REFERENCE); PRMs.emplace_back(PRM_STRING, 1, &ReferenceName, &ReferenceNameDefault);

GAS_ReadRefImage::GAS_ReadRefImage(const SIM_DataFactory* factory) : GAS_SubSolver(factory)
{
}

void GAS_ReadRefImage::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
}

void GAS_ReadRefImage::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
}

const SIM_DopDescription* GAS_ReadRefImage::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_REFERENCE

    static PRM_Name RefImageName("RefImage", "Reference Image");
    static PRM_Default RefImageDefault(0, "");
    PRMs.emplace_back(PRM_PICFILE, 1, &RefImageName, &RefImageDefault);

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

bool GAS_ReadRefImage::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_ScalarField* R = getScalarField(obj, GAS_NAME_REFERENCE);

    if (!R)
    {
        addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
        return false;
    }

    auto image = getRefImage();

    return true;
}
