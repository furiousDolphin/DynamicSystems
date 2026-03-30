
#ifndef SYSTEM_HPP_
#define SYSTEM_HPP_

#include <boost/math/tools/roots.hpp>
#include <boost/cstdint.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Dense>
#include <complex>

#include "ChebFunc.hpp"

#include "Settings.hpp"
#include "ValueManager.hpp"

class System
{
    public:
        System(Eigen::VectorXd b_coeffs, Eigen::VectorXd a_coeffs);

        using FuncType = std::variant<std::function<double(void)>, std::function<double(double)>>;

        std::array<Eigen::VectorXd, 2> step_response(double total_time, int n);
        std::array<Eigen::VectorXd, 2> impulse_response(double total_time, int n);

        std::array<Eigen::VectorXd, 3> nyquist_data(double D1, double D2, int n);
        std::array<Eigen::VectorXd, 3> bode_data(double D1, double D2, int n);


        void set_forcing_func(FuncType func);
        std::pair<double, double> do_RK4_step(double dt = 0.001);
        void operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double t);

    protected:
        void virtual update();

        enum class ResponseType
        { STEP, IMPULSE };
        std::array<Eigen::VectorXd, 2> step_impulse_helper(double total_time, int n, ResponseType response_type);

        enum class CharacteristicType
        { NYQUIST, BODE };
        std::array<Eigen::VectorXd, 3> nyquist_bode_helper(double D1, double D2, int n, CharacteristicType characteristic_type);


        struct TransferFunc
        {
            Eigen::VectorXd b_coeffs;
            Eigen::VectorXd a_coeffs;
            Eigen::MatrixXd M;
            void create_M();
        } tf_;

        struct ForcingFunc
        {
            std::function<double(double)> func;
            inline double operator()(double t);
        } forcing_func_;

        
        boost::numeric::odeint::runge_kutta4<Eigen::VectorXd> stepper_;
        Eigen::VectorXd state_;
        double t_;
        double dt_;
};


class SecondOrderSystem : public System
{
    public:
        SecondOrderSystem();

        using DoubleSetter = std::function<void(double)>;
        using DoubleGetter = std::function<double(void)>;

        struct Params
        {
            ValueManager zeta;
            ValueManager r;
            ValueManager f;
        };

        const Params& get_params() const;

    private:
        void update() override;
        Params params_;
};

#endif