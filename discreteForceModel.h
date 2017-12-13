#ifndef DISCRETEFORCEMODEL_H
#define DISCRETEFORCEMODEL_H

#include "spacecraft.h"
#include "integrationSettings.h"

class DiscreteForceModel
{
public:
    DiscreteForceModel(
            SpacecraftPointer spacecraftPointer,
            IntegrationSettingsPointer integrationSettingsPointer,
            ConstantsPointer constantsPointer);

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );

    ~DiscreteForceModel();

private:
    SpacecraftPointer spacecraftPointer_;
    IntegrationSettingsPointer integrationSettingsPointer_;
    ConstantsPointer constantsPointer_;
};

typedef boost::shared_ptr< DiscreteForceModel > DiscreteForceModelPointer;

#endif // DISCRETEFORCEMODEL_H
