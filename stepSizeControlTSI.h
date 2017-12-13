#ifndef STEPSIZECONTROLTSI_H
#define STEPSIZECONTROLTSI_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "integrationSettings.h"

class StepSizeControlTSI
{
public:
    StepSizeControlTSI( IntegrationSettingsPointer integrationSettingsPointer );

    ~StepSizeControlTSI();

    double determineNextStepSize( Eigen::MatrixXd derivatives, double previousStepSize );

private:
    IntegrationSettingsPointer integrationSettingsPointer_;

};

typedef boost::shared_ptr< StepSizeControlTSI > StepSizeControlTSIPointer;

#endif // STEPSIZECONTROLTSI_H
