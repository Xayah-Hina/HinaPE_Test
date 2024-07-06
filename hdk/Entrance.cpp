#include <UT/UT_DSOVersion.h> // Very Important!!! Include this first

#include "GAS_TestCubbyFlowSmoke.h"
#include "GAS_TestCubbyFlowFLIP.h"
#include "GAS_TestPhiFlow.h"

void initializeSIM(void *)
{
	IMPLEMENT_DATAFACTORY(GAS_TestCubbyFlowSmoke)
	IMPLEMENT_DATAFACTORY(GAS_TestCubbyFlowFLIP)
	IMPLEMENT_DATAFACTORY(GAS_TestPhiFlow)
}
