#include "ica.h"

//Approimating Functions
double gCosh(double x){
    return tanh(x);
}
double g_derCosh(double x){
    return 1.0 - std::pow(tanh(x),2.0);
}
double gGauss(double x){
    return x*exp(x*x/-2.0);
}
double g_derGauss(double x){
    return (1.0-x*x)*exp(x*x/-2.0);
}
double gCubic(double x){
    return .25*pow(x,4);;
}
double g_derCubic(double x){
    return pow(x,3);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructors
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set source Matrix
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
    Mat XX = SS*A.transpose();
    setSource(XX.transpose());
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Set initial mixing matrix
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Whiten data
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

//Gram-Schmidt Decorrelation for serial deflationary fastICA
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

//Symmetrical Decorrelation for parallel fastICA
void Ica::symDecorr(Mat* w)
{
    Eigen::EigenSolver<Mat> es(*w * w->transpose());
    Eigen::VectorXd D = es.eigenvalues().real();
    Mat V = es.eigenvectors().real();

    Mat D_to_neg_half = D.cwiseInverse().cwiseSqrt().asDiagonal();
    *w = V*D_to_neg_half*V.transpose() * *w;
}

//Based on https://www.cs.helsinki.fi/u/ahyvarin/papers/TNN99new.pdf
void Ica::serialFastICA(int n_sigs, double tol, int max_iter, std::function<double(double)> g, std::function<double(double)> g_der)
{
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

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

//Based on https://journal.r-project.org/archive/2018/RJ-2018-046/RJ-2018-046.pdf
void Ica::parallelFastICA(int n_sigs, double tol, int max_iter, std::function<double(double)> g, std::function<double(double)> g_der)
//Best seed 5
{
    std::cout << "In parallel" << std::endl;
    int M = X->cols(); 
    double M_inv = 1./M;
    Mat col_m = Mat::Ones(M,1); 

    symDecorr(W.get());

    for(int i=0; i<max_iter; i++)
    {
        Mat W_p = W->transpose() * *X; //W^T*X
        Mat E1 = M_inv * *X * W_p.unaryExpr(g).transpose(); //E{X*g(W^T*X)^T} 
        Mat E2 = *W; 
        Mat E2b = M_inv*(W_p.unaryExpr(g_der)*col_m); //E{g'(W^T*X)} 
        std::for_each(E2.colwise().begin(),E2.colwise().end(), //W.row(i) = E2b(i) * W.row(i)
            [&](auto&& col){            
                col.array() *= E2b.array(); 
            }
        );
        
        //Temp W
        auto Wa = std::make_unique<Mat>(E1 - E2);
        symDecorr(Wa.get());

        //Correlation between corresponding vectors in W and Wa, subtract 1 and get maximum correlation error
        double max_err = ((Wa->transpose() * *W).diagonal().cwiseAbs() - Mat::Ones(n_sigs,1)).cwiseAbs().maxCoeff(); 
        
        W = std::move(Wa);
        if (max_err < tol)
            break;

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

    //Func ptrs
    std::function<double(double)> g;
    std::function<double(double)> g_der;

    //Set func ptr to relavent approximating function
    if(func_type=="cosh"){
        g = gCosh;
        g_der = g_derCosh;
    }else if (func_type=="exp"){
        g = gGauss;
        g_der = g_derGauss;
    }else if (func_type=="cubic"){
        g = gCubic;
        g_der = g_derCubic;
    }else
        throw std::invalid_argument(func_type + " is not a valid function type option");

    //Call parallel or serial version
    if (paral == true)
        parallelFastICA(n_sigs, tol, max_iter,g,g_der);
    else{
        serialFastICA(n_sigs, tol, max_iter, g, g_der);
    }
}


