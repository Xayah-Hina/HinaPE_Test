#include "cubby.h"

void EmitSourcePartial(std::shared_ptr<CubbyFlow::ScalarGrid3>& TARGET, const SIM_RawField* SOURCE, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(SOURCE->field());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        float old_value = TARGET->DataView()(vit.x(), vit.y(), vit.z());
        TARGET->DataView()(vit.x(), vit.y(), vit.z()) = std::max(old_value, vit.getValue());
    }
}

THREADED_METHOD2(, true, EmitSource, std::shared_ptr<CubbyFlow::ScalarGrid3> &, TARGET, const SIM_RawField *, SOURCE);
void CubbyFlow::EmitFromHoudiniSource(std::shared_ptr<ScalarGrid3>& TARGET, const SIM_ScalarField* SOURCE) { EmitSource(TARGET, SOURCE->getField()); }


void EmitSourcePartial(std::shared_ptr<CubbyFlow::FaceCenteredGrid3>& TARGET, const SIM_RawField* SOURCE, int AXIS, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(SOURCE->field());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        float old_value = TARGET->DataView(AXIS)(vit.x(), vit.y(), vit.z());
        TARGET->DataView(AXIS)(vit.x(), vit.y(), vit.z()) = std::max(old_value, vit.getValue());
    }
}

THREADED_METHOD3(, true, EmitSource, std::shared_ptr<CubbyFlow::FaceCenteredGrid3> &, TARGET, const SIM_RawField *, SOURCE, int, AXIS);
void CubbyFlow::EmitFromHoudiniSource(std::shared_ptr<FaceCenteredGrid3>& TARGET, const SIM_VectorField* SOURCE) { for (int axis : {0, 1, 2}) { EmitSource(TARGET, SOURCE->getField(axis), axis); } }

void WriteFieldPartial(SIM_RawField* TARGET, const std::shared_ptr<CubbyFlow::ScalarGrid3>& SOURCE, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setArray(TARGET->fieldNC());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        vit.setValue((*SOURCE).DataView()(vit.x(), vit.y(), vit.z()));
    }
}

THREADED_METHOD2(, true, WriteField, SIM_RawField *, TARGET, const std::shared_ptr<CubbyFlow::ScalarGrid3> &, SOURCE);
void CubbyFlow::WriteHoudiniField(SIM_ScalarField* TARGET, const std::shared_ptr<ScalarGrid3>& SOURCE) { WriteField(TARGET->getField(), SOURCE); }

void WriteFieldPartial(SIM_RawField* TARGET, int axis, const std::shared_ptr<CubbyFlow::FaceCenteredGrid3>& SOURCE, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setArray(TARGET->fieldNC());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        vit.setValue(SOURCE->DataView(axis)(vit.x(), vit.y(), vit.z()));
    }
}

THREADED_METHOD3(, true, WriteField, SIM_RawField *, TARGET, int, axis, const std::shared_ptr<CubbyFlow::FaceCenteredGrid3> &, SOURCE);
void CubbyFlow::WriteHoudiniField(SIM_VectorField* TARGET, const std::shared_ptr<FaceCenteredGrid3>& SOURCE) { for (int axis : {0, 1, 2}) { WriteField(TARGET->getField(axis), axis, SOURCE); } }
