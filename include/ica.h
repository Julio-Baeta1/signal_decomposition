#ifndef ICA_H
#define ICA_H

#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <functional>
#include <utility>

using Mat = Eigen::MatrixXd;

class Ica{
    private:
        std::shared_ptr<Mat> X;
        std::unique_ptr<Mat> W;

        void gsDecorr(Mat *w, int col_num);
        void serialFastICACosh(int n_sigs, double tol, int max_iter);
        void serialFastICAExp(int n_sigs, double tol, int max_iter);
        void serialFastICACubic(int n_sigs, double tol, int max_iter);

    public:

        Ica(std::shared_ptr<Mat> ini_X = nullptr);

        void setSource(std::shared_ptr<Mat> new_X);
        void setSourceFromFile(std::string filename);
        std::shared_ptr<Mat> getResultPtr() {return X;}

        void randW(int rows, int cols, int seed);
        void setW(int rows, int cols, int seed, bool is_rand);

        void sphering();
        void decompose(int n_sigs=2, bool rand_W=false, int seed=1);
        void fastIca(int n_sigs=2, std::string func_type="cosh", int seed=1, double tol=1e-6, int max_iter=200);
};

#endif