#include "LKH.h"

void SOP_Report(GainType Cost)
{
    printff("Instance Cost = " GainFormat "_" GainFormat "\n",
            CurrentPenalty, Cost);
}
