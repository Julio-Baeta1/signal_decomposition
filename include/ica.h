#ifndef ICA_H
#define ICA_H

#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

using Mat = Eigen::MatrixXd;

class Ica{
    private:
        std::shared_ptr<Mat> X;
        std::unique_ptr<Mat> W;
    public:

        Ica(std::shared_ptr<Mat> ini_X = nullptr);

        void setSource(std::shared_ptr<Mat> new_X);
        void setSourceFromFile(std::string filename);
        std::shared_ptr<Mat> getResultPtr() {return X;}

        void randW(int rows, int cols, int seed);
        void setW(int rows, int cols, int seed, bool is_rand);

        void sphering();
        void decompose(int n_sigs, bool rand_W, int seed);
        void fastIca(int n_sigs, std::string func_type, int seed);
};

#endif