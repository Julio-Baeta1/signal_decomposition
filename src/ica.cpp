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

void Ica::decompose(int n_sigs)
/* Gradient Descent of Entropy
*/
{
    double n_samps = 2/(double)S->cols();
    
    /*std::mt19937 r_gen (123);
    std::normal_distribution<double> r_dis(0.0, 1.0);
    auto norm = [&] () {return r_dis(r_gen);};

    Mat W = Mat::NullaryExpr(n_sigs,n_sigs, norm );*/

    Mat W = Mat::Identity(n_sigs,n_sigs);
    Mat new_W = Mat::Zero(n_sigs,n_sigs);
    Mat u = Mat::Zero(S->rows(),n_samps);
    Mat U = Mat::Zero(S->rows(),n_samps);
    int max_iter = 1000;
    double step = 0.25;

    for(int i=0; i<max_iter; i++)
    {
        u = W * *S;
        U = u.array().tanh();
        new_W = n_samps* U * S->transpose();
        new_W = W.transpose().inverse() - new_W; 
        W = W + step*new_W;
    }
    std::cout << W.rows() <<" , " << W.cols() << std::endl << W << std::endl;
    *S = W* *S;
}