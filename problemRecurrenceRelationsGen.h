#ifndef PROBLEMRECURRENCERELATIONSGEN_H
#define PROBLEMRECURRENCERELATIONSGEN_H

#include "Eigen/Core"

#include "spacecraft.h"
#include "constants.h"


class ProblemRecurrenceRelations
{
public:
    ProblemRecurrenceRelations( bool useCSI,
                                bool useImprovedMass,
                                bool useSeparateSum,
                                int order,
                                SpacecraftPointer spacecraftPointer,
                                ConstantsPointer constantsPointer );

    ~ProblemRecurrenceRelations();

    int getOrder() {return order_;}

    void computeNextState(double currentStepSize);
//    Eigen::MatrixXd getNextStateVector() {return nextStateVector_;}
//    Eigen::MatrixXd getCoefficientMatrix() {return coefficientMatrix_;}

//    void setNextStateVector( Eigen::MatrixXd nextStateVector ) {nextStateVector_ = nextStateVector;}
//    void setCoefficientMatrix( Eigen::MatrixXd coefficientMatrix ) {coefficientMatrix_ = coefficientMatrix;}

private:
    bool useCSI_;
    bool useImprovedMass_;
    bool useSeparateSum_;
    int order_;
    SpacecraftPointer spacecraftPointer_;
    ConstantsPointer constantsPointer_;

//    Eigen::MatrixXd nextStateVector_;
//    Eigen::MatrixXd coefficientMatrix_;
};

typedef boost::shared_ptr< ProblemRecurrenceRelations > ProblemRecurrenceRelationsPointer;

#endif // PROBLEMRECURRENCERELATIONSGEN_H
