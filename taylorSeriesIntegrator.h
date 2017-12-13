#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <Eigen/Core>
#include "Eigen/Core"
#include "integrationSettings.h"
#include "spacecraft.h"
#include "recurrenceRelations.h"
#include "problemRecurrenceRelations.h"
#include "stepSizeControlTSI.h"

class TaylorSeriesIntegrator
{
public:
    TaylorSeriesIntegrator(
            ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer,
            SpacecraftPointer spacecraftpointer,
            IntegrationSettingsPointer integrationSettingsPointer,
            ConstantsPointer constantsPointer,
            StepSizeControlTSIPointer stepSizeControlTSIPointer);

    ~TaylorSeriesIntegrator();

    void integrate();

    void setStateMatrix( Eigen::MatrixXd stateMatrix ) {stateMatrix_ = stateMatrix;}
    void setTimeVector( Eigen::MatrixXd timeVector ) {timeVector_ = timeVector;}
    void setComputedThrustAccelerationMatrix( Eigen::MatrixXd computedThrustAccelerationMatrix )
            {computedThrustAccelerationMatrix_ = computedThrustAccelerationMatrix;}

    Eigen::MatrixXd getStateMatrix() {return stateMatrix_;}
    Eigen::MatrixXd getTimeVector() {return timeVector_;}
    Eigen::MatrixXd getComputedThrustAccelerationMatrix() {return computedThrustAccelerationMatrix_;}

/// DE VARIABELEN DIE DOOR DE GECOMMENTE FUNCTIE HIERONDER WORDT GEBRUIKT, ZIJN DEEL VAN DE SETTINGS OF
/// SETTINGS CLASS, DAARDOOR MOET IK ZE NIET ALLEMAAL APART MEEGEVEN MAAR KAN IK GEWOON
/// settingsPointer->maximumStepsize of spacecraftPointer->getCurrentAcceleration
//    void integratorTSI( char functionToEvaluate, double initialTime, double finalTime,
//                        Eigen::Matrix< double, 11, 1 > initialState, int order,
//                        int numberOfVariables,
//                        double minimumStepSize, double maximumStepSize,
//                        double initialStepSize, int refineFactor,
//                        thrustForceMatrixTypeDef thrustForceMatrix );
    std::vector<Vector11d> getStates() {return states_;}

private:
    void integratorStep(double& stepSize,const Eigen::Vector3d& thrustForceVector,const double limitingError);
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer_;
    SpacecraftPointer spacecraftPointer_;
    IntegrationSettingsPointer integrationSettingsPointer_;
    ConstantsPointer constantsPointer_;
    StepSizeControlTSIPointer stepSizeControlTSIPointer_;

    Eigen::MatrixXd stateMatrix_;
    Eigen::MatrixXd timeVector_;
    Eigen::MatrixXd computedThrustAccelerationMatrix_;

    std::vector<Vector11d> states_;


};

typedef boost::shared_ptr< TaylorSeriesIntegrator > TaylorSeriesIntegratorPointer;

#endif // INTEGRATOR_H
