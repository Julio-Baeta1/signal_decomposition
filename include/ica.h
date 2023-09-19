#ifndef ICA_H
#define ICA_H

#include <memory>
#include<Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/SVD>

using Mat = Eigen::MatrixXd;

class Ica{
    private:
        std::shared_ptr<Mat> S;

        void whiten();
    public:

        Ica(std::shared_ptr<Mat> initial_S = nullptr);

        void setSource(std::shared_ptr<Mat> new_S);
        std::shared_ptr<Mat> getResultPtr() {return S;}

        void decompose();
};

#endif