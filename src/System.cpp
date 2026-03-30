
#include <numbers>
#include <complex>
#include <unsupported/Eigen/Polynomials>

#include "System.hpp"



System::System(
    Eigen::VectorXd b_coeffs, 
    Eigen::VectorXd a_coeffs
) : 
    tf_{b_coeffs, a_coeffs},
    t_{0.0},
    dt_{0.0002}
{

}

void System::TransferFunc::create_M()
{
    int m = b_coeffs.size()-1;
    int n = a_coeffs.size()-1;
    if ( m > n )
    { throw std::runtime_error("m nie moze byc wieksze od n"); }

    M = Eigen::MatrixXd::Zero(n, n);
    if (n > 1)
    { M.block(0, 1, n-1, n-1) = Eigen::MatrixXd::Identity(n-1, n-1); }
    M.row(n-1) = a_coeffs.segment(0, n)*(-1/a_coeffs[n]);    
}

void System::update()
{
    tf_.create_M();
}

std::array<Eigen::VectorXd, 2> System::step_impulse_helper(double total_time, int n, ResponseType response_type)
{
    this->update();
    double a_0 = tf_.a_coeffs[0];
    int eq_order = tf_.a_coeffs.size()-1;
    state_ = Eigen::VectorXd::Zero(eq_order);

    switch (response_type)
    {
        case ResponseType::IMPULSE:
            forcing_func_.func = [](double t){return 0.0;};
            state_[eq_order-1] = 1/a_0;
            break;
        case ResponseType::STEP:
            forcing_func_.func = [](double t){return 1.0;};
            break;       
    }

    Eigen::VectorXd t_dense = Eigen::VectorXd::LinSpaced(n, 0.0, total_time);
    Eigen::VectorXd y_dense = Eigen::VectorXd::Zero(n);

    double dt = total_time/n;
    for (int idx = 0; idx<n; idx++)
    {
        auto [cur_t, cur_y] = this->do_RK4_step(dt);
        y_dense[idx] = cur_y;
    }
    return {t_dense, y_dense};
}

std::array<Eigen::VectorXd, 3> System::nyquist_bode_helper(double D1, double D2, int n, CharacteristicType characteristic_type)
{
    this->update();
    double PI = std::numbers::pi;

    std::function<double(std::complex<double>)> f1 = nullptr;
    std::function<double(std::complex<double>)> f2 = nullptr;

    switch (characteristic_type)
    {
        case CharacteristicType::NYQUIST:
            f1 = [&](std::complex<double> z){return z.real();};
            f2 = [&](std::complex<double> z){return z.imag();};
            break;
        case CharacteristicType::BODE:
            f1 = [&](std::complex<double> z){return 20.0*std::log10(std::abs(z));};
            f2 = [&](std::complex<double> z){return (180.0/PI)*std::arg(z);};
            break;       
    }

    Eigen::ArrayXd powers = Eigen::ArrayXd::LinSpaced(n, D1, D2);
    double ln10 = Eigen::numext::log(10.0);

    Eigen::VectorXd w_dense = Eigen::exp(powers*ln10).matrix();
    Eigen::VectorXd series1 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd series2 = Eigen::VectorXd::Zero(n);

    std::complex<double> imag_unit{0.0, 1.0};

    for (int idx = 0; idx < n; idx++)
    {
        double cur_w = w_dense[idx];
        std::complex<double> numerator = Eigen::poly_eval(tf_.b_coeffs, imag_unit*cur_w);
        std::complex<double> denominator = Eigen::poly_eval(tf_.a_coeffs, imag_unit*cur_w);

        std::complex<double> val = numerator/denominator;
        series1[idx] = f1(val);
        series2[idx] = f2(val);
    }

    return {w_dense, series1, series2};
}

std::array<Eigen::VectorXd, 2> System::impulse_response(double total_time, int n)
{ return this->step_impulse_helper(total_time, n, ResponseType::IMPULSE); }

std::array<Eigen::VectorXd, 2> System::step_response(double total_time, int n)
{ return this->step_impulse_helper(total_time, n, ResponseType::STEP); }

std::array<Eigen::VectorXd, 3> System::nyquist_data(double D1, double D2, int n)
{ return this->nyquist_bode_helper(D1, D2, n, CharacteristicType::NYQUIST); }   

std::array<Eigen::VectorXd, 3> System::bode_data(double D1, double D2, int n)
{ return this->nyquist_bode_helper(D1, D2, n, CharacteristicType::BODE); }


void System::set_forcing_func(FuncType func)
{
    std::visit(overloaded{
        [&](std::function<double(void)> void_func)
        {forcing_func_.func = [void_func](double t){return void_func();};},
        [&](std::function<double(double)> double_func)
        {forcing_func_.func = double_func;}
    }, func);

    t_ = 0.0;
    this->update();
    int N = 16;
    int m = tf_.b_coeffs.size() - 1;
    int n = tf_.a_coeffs.size() - 1;

    Eigen::VectorXd X_derivatives_vect = Eigen::VectorXd::Zero(n + 1); 

    std::visit(overloaded{
        [&](std::function<double(void)> void_func)
        { 
            if (m < N-1) 
            {
                double L = 0.001;
                std::pair<double, double> cheb_range{0.0, L};
                Eigen::VectorXd cheb_nodes = ChebNodes(N, cheb_range.first, cheb_range.second);
                Eigen::VectorXd cheb_vals = cheb_nodes.unaryExpr([&](double t){return void_func();});
                Eigen::VectorXd cheb_coeffs = ChebCoeffs(cheb_vals);
                Eigen::MatrixXd D = DifferentiationOperator(N); 
                Eigen::VectorXd kth_deriv_coeffs = cheb_coeffs;

                for (int i = 0; i < n; i++) 
                {
                    X_derivatives_vect[i] = Clenshaw(kth_deriv_coeffs, cheb_range.first);
                    kth_deriv_coeffs = D * kth_deriv_coeffs; 
                }
            } 
            else 
            { throw std::runtime_error("..."); }
        },
        [&](std::function<double(double)> double_func)
        { X_derivatives_vect[0] = double_func(0.0); }
    }, func);

    Eigen::MatrixXd X_derivatives_matrix = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) 
    {
        int len = n - i;
        X_derivatives_matrix.col(i).segment(i, len) = X_derivatives_vect.segment(0, len);
    }

    Eigen::VectorXd Y_derivatives_vect = Eigen::VectorXd::Zero(n);
    if (n > 1) 
    { Y_derivatives_vect[0] = (tf_.b_coeffs[0]/tf_.a_coeffs[0]) * X_derivatives_vect[1]; }

    Eigen::MatrixXd lin_eq_sys_matrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
    b.segment(0, m+1) = tf_.b_coeffs;

    lin_eq_sys_matrix.row(0) = b.transpose() * tf_.M;
    for (int i = 1; i < n; i++) 
    { lin_eq_sys_matrix.row(i) = lin_eq_sys_matrix.row(i-1) * tf_.M; }

    Eigen::VectorXd out_vect = lin_eq_sys_matrix.col(n-1).transpose() * X_derivatives_matrix;

    out_vect *= (-1.0 / tf_.a_coeffs[n]);
    out_vect += Y_derivatives_vect;

    state_ = lin_eq_sys_matrix.colPivHouseholderQr().solve(out_vect);

}

std::pair<double, double> System::do_RK4_step(double dt)
{
    this->update();
    double remaining_time = dt;

    while (remaining_time >= dt_)
    {
        stepper_.do_step(*this, state_, t_, dt_);
        t_ += dt_;
        remaining_time -= dt_;
    }

    if (remaining_time > 0.0)
    {
        stepper_.do_step(*this, state_, t_, remaining_time);
        t_ += remaining_time;
    }

    double y = 0.0;
    for(int i = 0; i < tf_.a_coeffs.size()-1; ++i) 
    { y += tf_.b_coeffs[i] * state_[i]; }

    return {t_, y};
}

void System::operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double t)
{
    int n = tf_.a_coeffs.size()-1;
    double a_n = tf_.a_coeffs[n];
    double x = forcing_func_(t);

    dZdt = tf_.M*Z;
    dZdt[n-1] += x/a_n;
}

inline double System::ForcingFunc::operator()(double t)
{ return func(t); }


SecondOrderSystem::SecondOrderSystem():
    System(
        (Eigen::VectorXd{2} << 1.0, 0.0).finished(),
        (Eigen::VectorXd{3} << 1.0, 2.0, 1.0).finished())
{

}

void SecondOrderSystem::update()
{
    auto& [zeta_vm, r_vm, f_vm] = params_;

    if 
    ( 
        zeta_vm.check_and_reset_dirty() || 
        r_vm.check_and_reset_dirty() || 
        f_vm.check_and_reset_dirty() 
    )

    {
        double zeta = zeta_vm.get_val();
        double r = r_vm.get_val();
        double f = f_vm.get_val();

        double PI = std::numbers::pi;
        double w = 2*PI*f;

        double k1 = (2*zeta)/w;
        double k2 = 1/(w*w);
        double k3 = (zeta*r)/w;    

        tf_.b_coeffs[0] = 1.0;
        tf_.b_coeffs[1] = k3;

        tf_.a_coeffs[0] = 1.0;
        tf_.a_coeffs[1] = k1; 
        tf_.a_coeffs[2] = k2;

        tf_.create_M();
    }
}

const SecondOrderSystem::Params& SecondOrderSystem::get_params() const
{ return params_; }


