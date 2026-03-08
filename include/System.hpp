
#ifndef SYSTEM_HPP_
#define SYSTEM_HPP_

#include <boost/math/tools/roots.hpp>
#include <boost/cstdint.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Dense>

class System
{
    public:
        System(Eigen::VectorXd b_coeffs, Eigen::VectorXd a_coeffs);
        enum ForcingFunctionType
        {
            STEP,
            IMPULSE
        };
        void set_forcing_function(ForcingFunctionType forcing_function_type);
        void set_forcing_function(const double* new_x_ptr_);
        void set_forcing_function(std::function<double(double)> forcing_function);

        void do_step();

    private:
        void operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double dt);
        struct TransferFunction
        {
            Eigen::VectorXd b_coeffs;
            Eigen::VectorXd a_coeffs;
        } tf_;
        const double* x_ptr_;
};

#endif