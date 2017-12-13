#ifndef RECURRENCERELATIONS
#define RECURRENCERELATIONS

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

class RecurrenceRelations
{
public:
    RecurrenceRelations( int orderInLoopOfTaylorSeries ):
        orderInLoopOfTaylorSeries_( orderInLoopOfTaylorSeries ){}

    double multiplyRecursive(
            const Eigen::VectorXd &X_left,      // dim: (order+1) x 1
            const Eigen::VectorXd &U_left,
            const Eigen::VectorXd &X_right,
            const Eigen::VectorXd &U_right );

    double V1overV2(
            const Eigen::VectorXd &W,          // dim: (order+1) x 1
            const Eigen::VectorXd &V1,
            const Eigen::VectorXd &V2 );

    double divideRecursive(
            const Eigen::VectorXd &W,          // dim: (order+1) x 1
            const Eigen::VectorXd &U_num,
            const Eigen::VectorXd &U_denom,
            const Eigen::VectorXd &X_denom );

    double Um_over_Xn(const Eigen::VectorXd &W,          // dim: (order+1) x 1
            const Eigen::VectorXd &U_num,
            const Eigen::VectorXd &U_denom,
            const Eigen::VectorXd &X_denom );

    double X1xV2_Recursive(
            const Eigen::VectorXd &X_L,          // dim: (order+1) x 1
            const Eigen::VectorXd &U_L,
            const Eigen::VectorXd &V_R );

    double V1xV2(
            const Eigen::VectorXd &V1,          // dim: (order+1) x 1
            const Eigen::VectorXd &V2 );

    double X1xU2overX2_Recursive(
            const Eigen::VectorXd &W,          // dim: (order+1) x 1
            const Eigen::VectorXd &X_L,
            const Eigen::VectorXd &U_L,
            const Eigen::VectorXd &X_R,
            const Eigen::VectorXd &U_R);

private:
    int orderInLoopOfTaylorSeries_;

};

typedef boost::shared_ptr< RecurrenceRelations > RecurrenceRelationsPointer;

#endif // RECURRENCERELATIONS

