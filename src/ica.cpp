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
    //W->normalize(); //Frobenius Norm
    W->rowwise().normalize();
    
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

void Ica::fastIca(int n_sigs, std::string func_type, int seed)
//Serial implementation, must make more parallel
//best seed currently 14 but returns signals scaled by -1
{

    //Lambda functions must fix scope problem
    /*if(func_type=="cosh"){
        auto g = [=] (double x) {return tanh(x);};
        auto g_der = [=] (double x) {return 1.0 - std::pow(tanh(x),2.0);};
    }else if (func_type=="exp")
    {
        auto g = [=] (double x) {return x*exp(x*x/-2.0);};
        auto g_der = [=] (double x) {return (1.0-x*x)*exp(x*x/-2.0);};
    }else if (func_type=="cubic")
    {
        auto g = [=] (double x) {return .25*pow(x,4);};
        auto g_der = [=] (double x) {return pow(x,3);};
    }else{
        throw std::invalid_argument(func_type + " is not a valid function type option");
    }*/

    //Lambda functions
    auto g = [=] (double x) {return tanh(x);};
    auto g_der = [=] (double x) {return 1.0 - std::pow(tanh(x),2.0);};

    //whiten data
    sphering();

    int N = X->rows();
    int M = X->cols(); 
    int max_iter = 1000;

    //Gen Random W
    randW(N, n_sigs, seed);

    //used to get Expected Value
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    for(int n=0; n< n_sigs; n++)
    {
        Mat w_p = W->row(n);

        for(int i=0; i<max_iter; i++)
        {
            w_p = w_p * *X; //projection

            Mat f1 = M_inv * *X * w_p.unaryExpr(g).transpose(); //E{X*g(w*X)^T}
            Mat f2 = w_p.unaryExpr(g_der)*col_m; //1x1 matrix to scale w
            Mat f3 = M_inv * f2(0,0) * W->row(n) ; //E{g'(w*X)}*w}

            w_p = f1.transpose() - f3;  

            //Gram Schimdt
            Mat sum = Mat::Zero(1,N);
            for (int k=0; k<n; k++)
                {
                    Mat temp = w_p * W->row(k).transpose();
                    sum += temp(0,0) * W->row(k);
                } 

            w_p = w_p - sum;
            w_p.normalize();
        }
        W->row(n) = w_p;
    }

    std::cout << *W << std::endl;
    *X = *W * *X;
}