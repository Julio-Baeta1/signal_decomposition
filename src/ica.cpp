#include "ica.h"

Ica::Ica(std::shared_ptr<Mat> initial_S)
{
    S = initial_S;
}

void Ica::setSource(std::shared_ptr<Mat> new_S)
{
    S = new_S;
}


void Ica::sphering()
{
    //centering
    Eigen::VectorXd S_m = S->rowwise().mean();
    *S = S->array().colwise()-S_m.array();

    //Auto-covariance matrix
    Mat S_cov = *S * S->transpose();

    //Eigen Decomposition
    Eigen::EigenSolver<Mat> es(S_cov);
    Eigen::VectorXd D = es.eigenvalues().real();
    Mat V = es.eigenvectors().real();

    //Whitening
    Mat D_to_neg_half = D.cwiseInverse().cwiseSqrt().asDiagonal();
    *S = V*D_to_neg_half*V.transpose() * *S;
}

void Ica::decompose()
{
    int x=1;
}