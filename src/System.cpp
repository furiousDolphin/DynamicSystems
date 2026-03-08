
#include "System.hpp"


System::System(
    Eigen::VectorXd b_coeffs, 
    Eigen::VectorXd a_coeffs
) : 
    tf_{b_coeffs, a_coeffs}
{

}

void System::set_forcing_function(System::ForcingFunctionType forcing_function_type)
{

}

void System::set_forcing_function(const double* new_x_ptr_)
{

}

void System::set_forcing_function(std::function<double(double)> forcing_function)
{

}

void System::do_step()
{

}

void System::operator()(const Eigen::VectorXd& Z, Eigen::VectorXd& dZdt, double dt)
{
    
}