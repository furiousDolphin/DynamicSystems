
#include "System.hpp"


System::System(
    Eigen::VectorXd b_coeffs, 
    Eigen::VectorXd a_coeffs
) : 
    tf_{b_coeffs, a_coeffs}
{
    int m = tf_.b_coeffs.size()-1;
    int n = tf_.a_coeffs.size()-1;
    if ( m > n )
    { throw std::runtime_error("m nie moze byc wieksze od n"); }

    // tf_.b = Eigen::VectorXd::Zero(m+1);
    // tf_.b.segment(0, m+1) = tf_.b_coeffs;

    tf_.M = Eigen::MatrixXd::Zero(n, n);
    tf_.M.block(1, 1, n-1, n-1) = Eigen::MatrixXd::Identity(n-1, n-1);
    tf_.M.row(n-1) = tf_.a_coeffs.segment(0, n)*(-1/tf_.a_coeffs[n]);
}


Eigen::VectorXd System::step_response(const Eigen::VectorXd& t)
{

}

Eigen::VectorXd System::impulse_response(const Eigen::VectorXd& t)
{

}

void System::set_forcing_func(FuncType func)
{
    int N = 16;
    int m = tf_.b_coeffs.size() - 1;
    if ( m < N-1 )
    {
        double L = 0.001;
        std::pair<double, double> cheb_range{0.0, L};

        Eigen::VectorXd cheb_nodes = ChebNodes(N, cheb_range.first, cheb_range.second);
        Eigen::VectorXd cheb_vals = cheb_nodes.unaryExpr(func);
        Eigen::VectorXd cheb_coeffs = ChebCoeffs(cheb_vals);

        Eigen::MatrixXd X_derivatives{m};
        Eigen::MatrixXd D = DifferentiationOperator(N); 

        for ( int i = 0; i < m; i++ )  
        { 
            X_derivatives[i] = Clenshaw(cheb_coeffs, cheb_range.first);
            cheb_coeffs = D*cheb_coeffs;
        }


    }
    else 
    { throw std::runtime_error("trzeba dodac mozliwosc rozniczkowania wiecej niz 16 razy"); }


}

std::pair<double, double> System::do_RK4_step()
{

}

void System::operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double t)
{

}

double System::ForcingFunc::get_x(double t)
{
    std::visit(overloaded{
        [&](std::function<double(void)> void_func)
        { return void_func(); },
        [&](std::function<double(double)> double_func)
        { return double_func(t); }
    }, func);
}
