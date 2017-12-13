#ifndef DISCRETEFORCEMODELGEN_H
#define DISCRETEFORCEMODELGEN_H

#include "spacecraft.h"
#include "integrationSettings.h"

class DiscreteForceModel
{
public:
    DiscreteForceModel(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer);

    virtual void updateCurrentForcesAndAccelerationsForTSI(
                    double currentTime ) = 0;

    virtual ~DiscreteForceModel();
protected:
    SpacecraftPointer spacecraftPointer_;
    IntegrationSettingsPointer integrationSettingsPointer_;
    ConstantsPointer constantsPointer_;
};

typedef boost::shared_ptr< DiscreteForceModel > DiscreteForceModelPointer;

class DiscreteForceModelL2: public DiscreteForceModel
{
public:
    DiscreteForceModelL2(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer):DiscreteForceModel( spacecraftPointer, integrationSettingsPointer, constantsPointer ){ }

    ~DiscreteForceModelL2(){ }

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );
};

class DiscreteForceModelL3: public DiscreteForceModel
{
public:
    DiscreteForceModelL3(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer):DiscreteForceModel( spacecraftPointer, integrationSettingsPointer, constantsPointer ){ }

    ~DiscreteForceModelL3(){ }

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );
};

class DiscreteForceModelCSI: public DiscreteForceModel
{
public:
    DiscreteForceModelCSI(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer):DiscreteForceModel( spacecraftPointer, integrationSettingsPointer, constantsPointer ){ }

    ~DiscreteForceModelCSI(){ }

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );
};

class DiscreteForceModelCSIcorrectGamma: public DiscreteForceModel
{
public:
    DiscreteForceModelCSIcorrectGamma(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer):DiscreteForceModel( spacecraftPointer, integrationSettingsPointer, constantsPointer ){ }

    ~DiscreteForceModelCSIcorrectGamma(){ }

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );
};

class DiscreteForceModelCSIimprovedMass: public DiscreteForceModel
{
public:
    DiscreteForceModelCSIimprovedMass(
    SpacecraftPointer spacecraftPointer,
    IntegrationSettingsPointer integrationSettingsPointer,
    ConstantsPointer constantsPointer):DiscreteForceModel( spacecraftPointer, integrationSettingsPointer, constantsPointer ){ }

    ~DiscreteForceModelCSIimprovedMass(){ }

    void updateCurrentForcesAndAccelerationsForTSI(
            double currentTime );
};

#endif // DISCRETEFORCEMODELGEN_H
