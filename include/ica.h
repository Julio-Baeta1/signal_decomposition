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
    public:

        Ica(std::shared_ptr<Mat> initial_X = nullptr);

        void setSource(std::shared_ptr<Mat> new_X);
        void setSourceFromFile(std::string filename);
        std::shared_ptr<Mat> getResultPtr() {return X;}

        void sphering();
        void decompose(int n_sigs);
};

#endif