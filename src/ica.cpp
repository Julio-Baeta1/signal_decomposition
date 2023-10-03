#include "ica.h"

Ica::Ica(std::shared_ptr<Mat> ini_X)
{
    X = std::move(ini_X);
    W = nullptr;
}

void Ica::setSource(std::shared_ptr<Mat> new_X)
{
    X = std::move(new_X);
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
    std::shared_ptr<Mat> XX = std::make_shared<Mat>((SS*A.transpose()).transpose());
    setSource(XX);
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
        u = *W * *X;
        U = u.array().tanh();
        new_W = norm_const* *X * U.transpose();
        new_W = W->transpose().inverse() - new_W; 
        *W = *W + step*new_W;
    }
    *X = *W * *X;     
}

void Ica::fastIca(int n_sigs)
//Serial implementation, must make more parallel
{
    sphering();

    int N = X->rows();
    int M = X->cols();
    double M_inv = 1./M; 
    int max_iter = 10000;

    randW(N, n_sigs, 3);
    Mat col_m = Mat::Ones(M,1); //col vec

    auto g = [=] (double x) {return tanh(x);};
    auto g_der = [=] (double x) {return 1.0 - std::pow(tanh(x),2.0);};

    for(int n=0; n< n_sigs; n++)
    {
        for(int i=0; i<max_iter; i++)
        {
            Mat w_p = W->col(n).transpose() * *X;

            Mat f1 = M_inv * *X * w_p.unaryExpr(g).transpose();
            Mat f2 = w_p.unaryExpr(g_der)*col_m;
            double f2_d = f2(0,0);
            Mat f3 = M_inv * f2_d * W->col(n) ;

            w_p = f1 - f3;

            Mat sum = Mat::Zero(N,1);
            for (int k=0; k<n; k++)
                {
                    Mat temp = w_p.transpose() * W->col(k);
                    double temp_d = temp(0,0);
                    sum += temp_d * W->col(k);
                } 

            w_p = w_p - sum;
            w_p.normalize();
            W->col(n) = w_p;
        }
    }
    *X = *W * *X;
}