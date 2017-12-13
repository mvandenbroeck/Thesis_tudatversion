#ifndef TAYLORSERIESINTEGRATORGENMAP_H
#define TAYLORSERIESINTEGRATORGENMAP_H

#include <vector>
#include <map>
#include <Eigen/Core>
#include "Eigen/Core"
#include "integrationSettings.h"
#include "spacecraft.h"
#include "recurrenceRelations.h"
#include "problemRecurrenceRelationsGen.h"
#include "discreteForceModelGen.h"
#include "stepSizeControlTSI.h"

class TaylorSeriesIntegrator
{
public:
    TaylorSeriesIntegrator(
            bool useCSI,
            bool useCorrectGamma,
            ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer,
            DiscreteForceModelPointer discreteForceModelPointer,
            SpacecraftPointer spacecraftPointer,
            IntegrationSettingsPointer integrationSettingsPointer,
            ConstantsPointer constantsPointer,
            StepSizeControlTSIPointer stepSizeControlTSIPointer);

    ~TaylorSeriesIntegrator( );

    void integrate();

    void setStateMap( std::map< double, Eigen::MatrixXd > stateMap ) {stateMap_ = stateMap;}
    void setComputedThrustAccelerationMap( std::map< double, Eigen::MatrixXd > computedThrustAccelerationMap )
            {computedThrustAccelerationMap_ = computedThrustAccelerationMap;}
    void setFullStateMap( std::map< double, Eigen::MatrixXd > fullStateMap )
            {fullStateMap_ = fullStateMap;}

   std::map< double, Eigen::MatrixXd > getStateMap() {return stateMap_;}
   std::map< double, Eigen::MatrixXd > getComputedThrustAccelerationMap() {return computedThrustAccelerationMap_;}
   std::map< double, Eigen::MatrixXd > getFullStateMap() {return fullStateMap_;}

   std::vector<Vector11d> getStates() {return states_;}

private:
    bool useCSI_;
    bool useCorrectGamma_;
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer_;
    DiscreteForceModelPointer discreteForceModelPointer_;
    SpacecraftPointer spacecraftPointer_;
    IntegrationSettingsPointer integrationSettingsPointer_;
    ConstantsPointer constantsPointer_;
    StepSizeControlTSIPointer stepSizeControlTSIPointer_;

    std::map< double, Eigen::MatrixXd > stateMap_;
    std::map< double, Eigen::MatrixXd > computedThrustAccelerationMap_;
    std::map< double, Eigen::MatrixXd > fullStateMap_;

    std::vector<Vector11d> states_;


};

typedef boost::shared_ptr< TaylorSeriesIntegrator > TaylorSeriesIntegratorPointer;

#endif // TAYLORSERIESINTEGRATORGEN_H
