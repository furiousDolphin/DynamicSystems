
#include "System.hpp"


System::System(
    Eigen::VectorXd b_coeffs, 
    Eigen::VectorXd a_coeffs
) : 
    tf_{b_coeffs, a_coeffs},
    t_{0.0},
    dt_{0.0001}
{
    int m = tf_.b_coeffs.size()-1;
    int n = tf_.a_coeffs.size()-1;
    if ( m > n )
    { throw std::runtime_error("m nie moze byc wieksze od n"); }

    tf_.M = Eigen::MatrixXd::Zero(n, n);
    tf_.M.block(1, 1, n-1, n-1) = Eigen::MatrixXd::Identity(n-1, n-1);
    tf_.M.row(n-1) = tf_.a_coeffs.segment(0, n)*(-1/tf_.a_coeffs[n]);
}


Eigen::VectorXd System::step_response(const Eigen::VectorXd& t)
{
    forcing_func_.func = [](){return 1.0;};

    int n = tf_.a_coeffs.size() - 1;
    double a_0 = tf_.a_coeffs[0];
    state_ = Eigen::VectorXd::Zero(n);
    state_[0] = 1.0/a_0;

    for (int i = 1; i < t.size(); i++)
    { this->do_RK4_step(t[i]-t[i-1]); }
}

Eigen::VectorXd System::impulse_response(const Eigen::VectorXd& t)
{
    int n = tf_.a_coeffs.size()-1;
    double a_n = tf_.a_coeffs[n];
    state_ = Eigen::VectorXd::Zero(n);
    state_[n-1] = 1.0/a_n;

    for (int i = 1; i < t.size(); i++)
    { this->do_RK4_step(t[i]-t[i-1]); }
}

void System::set_forcing_func(FuncType func)
{
    forcing_func_.func = func;

    if ( auto* func_ptr = std::get_if<std::function<double(double)>*>(func) )
    {

    }

    t_ = 0.0; //zeruje t

    int N = 16;
    int m = tf_.b_coeffs.size() - 1;
    int n = tf_.a_coeffs.size() - 1;

    if ( m < N-1 )
    {
        double L = 0.001;
        std::pair<double, double> cheb_range{0.0, L};

        Eigen::VectorXd cheb_nodes = ChebNodes(N, cheb_range.first, cheb_range.second);
        Eigen::VectorXd cheb_vals = cheb_nodes.unaryExpr(func);
        Eigen::VectorXd cheb_coeffs = ChebCoeffs(cheb_vals);

        Eigen::MatrixXd D = DifferentiationOperator(N); 
        Eigen::VectorXd kth_deriv_coeffs = cheb_coeffs;

        Eigen::VectorXd X_derivatives_vect(n); // wektor x(0), x'(0)...x^(n-1)(0)
        for (int i = 0; i < n; ++i) 
        {
            X_derivatives_vect[i] = Clenshaw(kth_deriv_coeffs, cheb_range.first);
            kth_deriv_coeffs = D * kth_deriv_coeffs; 
        }

        // Eigen::VectorXd X_derivatives_vect = Eigen::VectorXd::Zero(n-1);
        // X_derivatives_vect[0] = forcing_func_(0.0);       

        Eigen::MatrixXd X_derivatives_matrix = Eigen::MatrixXd::Zero(n, n);

        Eigen::VectorXd Y_derivatives_vect = Eigen::VectorXd::Zero(n);
        Y_derivatives_vect[0] = (tf_.b_coeffs[0]/tf_.a_coeffs[0]) * X_derivatives_vect[1];

        Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
        b.segment(0, m+1) = tf_.b_coeffs;

        Eigen::VectorXd out_vect = Eigen::VectorXd::Zero(n);
        Eigen::MatrixXd lin_eq_sys_matrix = Eigen::MatrixXd::Zero(n, n);

        lin_eq_sys_matrix.row(0) = tf_.M*b;
        for ( int i = 1; i < n; i++ )
        { lin_eq_sys_matrix.row(i) = tf_.M*lin_eq_sys_matrix.row(i-1); }

        for ( int i = 1; i < n; i++ )
        { X_derivatives_matrix.col(i).segment(i, n-i) = X_derivatives_vect.segment(0, n-i);}

        
        //out_vect = (X_derivatives_matrix * lin_eq_sys_matrix).col(n-1);  alternatywna wersja liczenia out_vect ( to jest tylko zamiast tego for )

        for (int i = 0; i < n; i++)
        {
            for ( int j = 0; j < n; j++ )
            { out_vect[j]+=lin_eq_sys_matrix.row(i)[n]*X_derivatives_matrix.row(i)[j]; }       
        }

        out_vect*=(-1.0/tf_.a_coeffs[n]);
        out_vect+=Y_derivatives_vect;

        state_ = lin_eq_sys_matrix.colPivHouseholderQr().solve(out_vect);
    }
    else 
    { throw std::runtime_error("trzeba dodac mozliwosc rozniczkowania wiecej niz 16 razy"); }


}

std::pair<double, double> System::do_RK4_step(double dt)
{
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
    for(int i = 0; i < tf_.b_coeffs.size(); ++i) 
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

double System::ForcingFunc::operator()(double t)
{
    return std::visit(overloaded{
        [&](std::function<double(void)> void_func)
        { return void_func(); },
        [&](std::function<double(double)> double_func)
        { return double_func(t); }
    }, func);
}
