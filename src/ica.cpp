#include "ica.h"

Ica::Ica(std::shared_ptr<Mat> initial_X)
{
    X = initial_X;
}

void Ica::setSource(std::shared_ptr<Mat> new_X)
{
    X = new_X;
}

void Ica::setSourceFromFile(std::string filename)
//Currently hardcoded to scikit-learn example
{
    Mat SS = Mat::Zero(2000,3);
    Mat A{{1, 1, 1}, {0.5, 2, 1.0}, {1.5, 1.0, 2.0}};

    std::ifstream mat_file (filename);
    std::string line, elem;
    int i{0};
    while (std::getline (mat_file, line)) {
        std::stringstream ss(line);
        int j{0};
        while (std::getline(ss, elem, ','))
        {
            SS(i,j) = stod(elem);
            j++;
        }
        i++;
    }

    *X = SS*A.transpose();
}


void Ica::sphering()
{
    //centering
    Eigen::VectorXd X_m = X->rowwise().mean();
    *X = X->array().colwise()-X_m.array();

    //Auto-covariance matrix
    Mat X_cov = *X * X->transpose();

    //Eigen Decomposition
    Eigen::EigenSolver<Mat> es(X_cov);
    Eigen::VectorXd D = es.eigenvalues().real();
    Mat V = es.eigenvectors().real();

    //Whitening
    Mat D_to_neg_half = D.cwiseInverse().cwiseSqrt().asDiagonal();
    *X = V*D_to_neg_half*V.transpose() * *X;
}

void Ica::decompose(int n_sigs)
/* Gradient Descent of Entropy
*/
{   
    //Mat X = S->transpose(); 

    int num_epoch = 10000;
    double learning_rate = 0.00003;

    Mat W2 = Mat::Identity(3,3);
    Mat new_W2 = Mat::Zero(3,3);
    Mat u2 = Mat::Zero(2000,3);
    Mat U2 = Mat::Zero(2000,3);

    for(int i=0; i < num_epoch; i++)
    {
        u2 = *X * W2;
        U2 = u2.array().tanh();
        new_W2 = W2.transpose().inverse() - (2.0/2000) * X->transpose() * U2;
        W2 = W2 + learning_rate*new_W2;
    }
    
    *X = *X * W2;
    //*S = est_S.transpose();
    
    //std::mt19937 r_gen (123);
    //std::normal_distribution<double> r_dis(0.0, 1.0);
    //auto norm = [&] () {return r_dis(r_gen);};

    //Mat W = Mat::NullaryExpr(n_sigs,n_sigs, norm );

    /*double n_samps = 2/(double)S->cols();
    Mat W = Mat::Identity(n_sigs,n_sigs);
    Mat new_W = Mat::Zero(n_sigs,n_sigs);
    Mat u = Mat::Zero(S->rows(),n_samps);
    Mat U = Mat::Zero(S->rows(),n_samps);
    int max_iter = 10000;
    double step = 0.00003;

    for(int i=0; i<max_iter; i++)
    {
        u = W * *S;
        U = u.array().tanh();
        new_W = n_samps* *S * U.transpose();
        new_W = W.transpose().inverse() - new_W; 
        W = W + step*new_W;
    }
    *S = W* *S;*/
}