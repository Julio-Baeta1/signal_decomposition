#include "ica.h"

Ica::Ica(Mat ini_X)
{
    X = std::make_unique<Mat>();
    W = std::make_unique<Mat>();
    XX = std::make_unique<Mat>(ini_X);
}

Ica::Ica(Mat* ini_ptr)
{
    X = std::make_unique<Mat>();
    W = std::make_unique<Mat>();
    XX = std::make_unique<Mat>(*ini_ptr);
}

void Ica::setSource(Mat new_X)
{
    XX = std::make_unique<Mat>(new_X);
}

void Ica::setSource(Mat* new_ptr)
{
    XX = std::make_unique<Mat>(*new_ptr);
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
    //std::shared_ptr<Mat> XX = std::make_shared<Mat>((SS*A.transpose()).transpose());
    Mat XX = SS*A.transpose();
    setSource(XX.transpose());
    //setSource(XX);
}

void Ica::randW(int rows, int cols, int seed)
{
    std::mt19937 r_gen (seed);
    std::normal_distribution<double> r_dis(0.0, 1.0);
    auto norm = [&] () {return r_dis(r_gen);};

    W = std::make_unique<Mat>(Mat::NullaryExpr(rows,cols, norm )) ;
    
    W->rowwise().normalize(); //W->normalize() uses the Frobenius Norm
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
    *X = D_to_neg_half*V.transpose() * *X;
}

void Ica::gsDecorr(Mat* w, int col_num)
{
    //Gram Schimdt
    Mat temp;
    Mat sum = Mat::Zero(X->rows(),1);
    for (int k=0; k<col_num; k++)
    {
        temp = w->transpose() * W->col(k);
        sum += temp(0,0) * W->col(k);
    } 

    *w = *w - sum;
    w->normalize();

}

void Ica::serialFastICACosh(int n_sigs, double tol, int max_iter)
{
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    auto g = [=] (double x) {return tanh(x);};
    auto g_der = [=] (double x) {return 1.0 - std::pow(tanh(x),2.0);};

    for(int n=0; n< n_sigs; n++)
    {
        Mat w_p = W->col(n);

        for(int i=0; i<max_iter; i++)
        {
            std::cout << "i: " <<n << std::endl << std::endl;
            Mat w_proj = w_p.transpose() * *X; //w^T*X
            Mat E1 = M_inv * *X * w_proj.unaryExpr(g).transpose(); //E{X*g(w^T*X)^T}
            Mat E2 = M_inv*(w_proj.unaryExpr(g_der)*col_m)(0,0) * w_p ; //E{g'(w^T*X)}*w} 
            w_proj = E1 - E2; 

            gsDecorr(&w_proj, n);

            //err = the absolute difference between 1 and the correlation between w_p(t) and w_p(t-1) 
            double err = abs(abs((w_proj.transpose() * w_p)(0,0)) -1); 
            w_p = w_proj;
            if(err < tol)
                break;
            
        }
        W->col(n) = w_p;
    }

    std::cout << W->transpose() << std::endl;
    *X = W->transpose() * *X;
}

void Ica::serialFastICAExp(int n_sigs, double tol, int max_iter)
{
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    auto g = [=] (double x) {return x*exp(x*x/-2.0);};
    auto g_der = [=] (double x) {return (1.0-x*x)*exp(x*x/-2.0);};

    for(int n=0; n< n_sigs; n++)
    {
    
        Mat w_p = W->col(n);

        for(int i=0; i<max_iter; i++)
        {
            Mat w_proj = w_p.transpose() * *X; //w^T*X
            Mat E1 = M_inv * *X * w_proj.unaryExpr(g).transpose(); //E{X*g(w^T*X)^T}
            Mat E2 = M_inv*(w_proj.unaryExpr(g_der)*col_m)(0,0) * w_p ; //E{g'(w^T*X)}*w} 
            w_proj = E1 - E2; 

            gsDecorr(&w_proj, n);

            //err = the absolute difference between 1 and the correlation between w_p(t) and w_p(t-1) 
            double err = abs(abs((w_proj.transpose() * w_p)(0,0)) -1); 
            w_p = w_proj;
            if(err < tol)
                break;
            
        }
        W->col(n) = w_p;
    }

    std::cout << W->transpose() << std::endl;
    *X = W->transpose() * *X;
}

void Ica::serialFastICACubic(int n_sigs, double tol, int max_iter)
{
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    auto g = [=] (double x) {return .25*pow(x,4);};
    auto g_der = [=] (double x) {return pow(x,3);};

    for(int n=0; n< n_sigs; n++)
    {
    
        Mat w_p = W->col(n);

        for(int i=0; i<max_iter; i++)
        {
            Mat w_proj = w_p.transpose() * *X; //w^T*X
            Mat E1 = M_inv * *X * w_proj.unaryExpr(g).transpose(); //E{X*g(w^T*X)^T}
            Mat E2 = M_inv*(w_proj.unaryExpr(g_der)*col_m)(0,0) * w_p ; //E{g'(w^T*X)}*w} 
            w_proj = E1 - E2; 

            gsDecorr(&w_proj, n);

            //err = the absolute difference between 1 and the correlation between w_p(t) and w_p(t-1) 
            double err = abs(abs((w_proj.transpose() * w_p)(0,0)) -1); 
            w_p = w_proj;
            if(err < tol)
                break;
            
        }
        W->col(n) = w_p;
    }

    std::cout << W->transpose() << std::endl;
    *X = W->transpose() * *X;
}

void Ica::parallelFastICACosh(int n_sigs, double tol, int max_iter)
//Best seed 5
{
    std::cout << "In parallel" << std::endl;
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    Eigen::EigenSolver<Mat> es(*W * W->transpose());
    Eigen::VectorXd D = es.eigenvalues().real();
    Mat V = es.eigenvectors().real();

    Mat D_to_neg_half = D.cwiseInverse().cwiseSqrt().asDiagonal();
    *W = V*D_to_neg_half*V.transpose() * *W;

    auto g = [=] (double x) {return tanh(x);};
    auto g_der = [=] (double x) {return 1.0 - std::pow(tanh(x),2.0);};

    for(int i=0; i<10; i++)
    {
        for(int k=0; k < n_sigs; k++)
        {
            Mat w_p = W->col(k);
            Mat w_proj = w_p.transpose() * *X; //w^T*X
            Mat E1 = M_inv * *X * w_proj.unaryExpr(g).transpose(); //E{X*g(w^T*X)^T}
            Mat E2 = M_inv*(w_proj.unaryExpr(g_der)*col_m)(0,0) * w_p ; //E{g'(w^T*X)}*w} 
            W->col(k)= E1 - E2; 
        }

        Eigen::EigenSolver<Mat> es(*W *W->transpose());
        Eigen::VectorXd D = es.eigenvalues().real();
        Mat V = es.eigenvectors().real();

        Mat D_to_neg_half = D.cwiseInverse().cwiseSqrt().asDiagonal();
        *W = V*D_to_neg_half*V.transpose() * *W;
    }

    std::cout << W->transpose() << std::endl;
    *X = W->transpose() * *X;
}

void Ica::fastIca(int n_sigs, std::string func_type, bool paral, int seed, double tol, int max_iter)
{
    X = std::make_unique<Mat>(*XX);

    //whiten data
    sphering();

    //Gen Random W
    randW(X->rows(), n_sigs, seed);

    

    /*Call different sub functions with the appropriate lambda functions for the approximating function selected. A more elegant 
    solution would be to use function pointer std::function<double(double)> but comipler optimiser/linker issue caused some func ptr
    have undefined reference to the function. This will be resolved at a later stage as the current focus is on having a complete 
    working project prototype sooner rather than an elegant solution later.
    */

    if (paral == true)
        parallelFastICACosh(n_sigs, tol, max_iter);

    else if(func_type=="cosh")
        serialFastICACosh(n_sigs, tol, max_iter);
    else if (func_type=="exp")
        serialFastICAExp(n_sigs, tol, max_iter);
    else if (func_type=="cubic")
        serialFastICACubic(n_sigs, tol, max_iter);
    else
        throw std::invalid_argument(func_type + " is not a valid function type option");
}


