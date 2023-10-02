#include "ica.h"

Ica::Ica(std::shared_ptr<Mat> ini_X)
{
    X = std::move(ini_X);
    W = nullptr;//std::make_unique<Mat>(Mat::Identity(X->cols(),X->cols())) ;
}

void Ica::setSource(std::shared_ptr<Mat> new_X)
{
    std::cout << "Inside set Source" << std::endl;
    X = std::move(new_X);
}


void Ica::setSourceFromFile(std::string filename)
//Currently hardcoded to scikit-learn example
{
    std::cout << "Inside Source from file" << std::endl;
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
    std::shared_ptr<Mat> XX = std::make_shared<Mat>((SS*A.transpose()).transpose());
    setSource(XX);
    //X = std::make_shared<Mat>((SS*A.transpose()).transpose());
    //X = std::make_shared<Mat>(SS*A.transpose());
}

void Ica::randW(int rows, int cols, int seed)
{
    std::mt19937 r_gen (seed);
    std::normal_distribution<double> r_dis(0.0, 1.0);
    auto norm = [&] () {return r_dis(r_gen);};

    W = std::make_unique<Mat>(Mat::NullaryExpr(rows,cols, norm )) ;
}

void Ica::setW(int rows, int cols, int seed, bool is_rand)
{
    if(is_rand)
        randW(rows, cols, seed);
    else
        W = std::make_unique<Mat>(Mat::Identity(rows,cols)) ;
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

void Ica::decompose(int n_sigs, bool rand_W, int seed)
/* Gradient Descent of Entropy
*/
{   
    setW(n_sigs,n_sigs,seed,rand_W);

    double norm_const = 2/(double)X->cols();
    Mat new_W = Mat::Zero(n_sigs,n_sigs);
    Mat u = Mat::Zero(X->rows(),X->cols());
    Mat U = Mat::Zero(X->rows(),X->cols());
    int max_iter = 10000;
    double step = 0.00003;

    for(int i=0; i<max_iter; i++)
    {
        //u = W * *X;
        u = *W * *X;
        U = u.array().tanh();
        new_W = norm_const* *X * U.transpose();
        //new_W = W.transpose().inverse() - new_W; 
        //W = W + step*new_W;
        new_W = W->transpose().inverse() - new_W; 
        *W = *W + step*new_W;
    }
    *X = *W * *X;     
}

/*void Ica::fastIca(int n_sigs)
{
    sphering();
    Mat XX = X->transpose();

    int M = XX.rows(), N= XX.cols();
    int max_iter = 1000;

    Mat W

    for(int i=0; i<max_iter; i++)
    {

    }

}*/