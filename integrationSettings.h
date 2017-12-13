#ifndef INTEGRATIONSETTINGS_H
#define INTEGRATIONSETTINGS_H

#include "constants.h"
#include "Eigen/Core"

class IntegrationSettings
{
public:
    IntegrationSettings( double initialTime, double finalTime );
    double getInitialTime() {return initialTime_;}
    double getFinalTime() {return finalTime_;}
    double getSafetyFactor() {return safetyFactor_;}
    // Special function to acces one element in error tolerance at a time
    double getErrorTolerance( int index ) {return errorTolerance_( index );}
    double getInitialStepSize() {return initialStepSize_;}
    double getStepSize() {return stepSize_;}
    int getOrderOfTaylorSeries() {return orderOfTaylorSeries_;}
    double getMinimumStepSize() {return minimumStepSize_;}
    double getMaximumStepSize() {return maximumStepSize_;}
    int getRefineFactor() {return refineFactor_;}

    //void setInitialTime(double initialTime) {initialTime_ = initialTime;}
    //void setFinalTime(double finalTime) {finalTime_ = finalTime;}
    void setSafeftyFactor(double safetyFactor) {safetyFactor_ = safetyFactor;}
    void setErrorTolerance( Eigen::Matrix< double, Eigen::Dynamic, 1 > errorTolerance )
                            {errorTolerance_ = errorTolerance;}
    void setInitialStepSize(double initialStepsize) {initialStepSize_ = initialStepsize;}
    void setStepSize(double stepSize) {stepSize_ = stepSize;}
    void setMinimumStepSize( double minimumStepSize ) {minimumStepSize_ = minimumStepSize;}
    void setMaximumStepSize( double maximumStepSize ) {maximumStepSize_ = maximumStepSize;}
    void setRefineFactor( int refineFactor ) {refineFactor_ = refineFactor;}
    void setOrderOfTaylorSeries( int orderOfTaylorSeries ) {orderOfTaylorSeries_ = orderOfTaylorSeries;}

private:
    double initialTime_;
    double finalTime_;
    double safetyFactor_;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > errorTolerance_;
    double initialStepSize_;
    double stepSize_; // This stepsize will be constantly updated by integratorTSI
    double minimumStepSize_;
    double maximumStepSize_;
    int refineFactor_;
    int orderOfTaylorSeries_;
};

typedef boost::shared_ptr< IntegrationSettings > IntegrationSettingsPointer;

#endif // INTEGRATIONSETTINGS_H
