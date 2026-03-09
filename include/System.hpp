
#ifndef SYSTEM_HPP_
#define SYSTEM_HPP_

#include <boost/math/tools/roots.hpp>
#include <boost/cstdint.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Dense>
#include "ChebFunc.hpp"

#include "Settings.hpp"

class System
{
    public:
        System(Eigen::VectorXd b_coeffs, Eigen::VectorXd a_coeffs);

        using FuncType = std::variant<std::function<double(void)>, std::function<double(double)>>;

        Eigen::VectorXd step_response(const Eigen::VectorXd& t);
        Eigen::VectorXd impulse_response(const Eigen::VectorXd& t);

        void set_forcing_func(FuncType func);
        std::pair<double, double> do_RK4_step();

    private:
        void operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double t);

        struct TransferFunc
        {
            Eigen::VectorXd b_coeffs;
            Eigen::VectorXd a_coeffs;
            Eigen::MatrixXd M;
        } tf_;

        struct ForcingFunc
        {
            FuncType func;
            double get_x(double t);
        } forcing_func_;
};

#endif