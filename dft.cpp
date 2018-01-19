// translated from DFT Python code as shown in 
// https://www.youtube.com/watch?v=sX-DNi_SX-Q
//
// Compile:
// g++ -g --std=c++11 -o dft dft.cpp  -Ieigen
//
// Execute:
// ./dft  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
//
// Plot using gnuplot:
// plot "plot_out.txt" using 1:2 with lines lw 3, "plot_in.txt" using 1:2:(0.03) with circles fill solid

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>

using std::cout;
using namespace Eigen;
using std::pow;
using std::exp;
using std::ofstream;
using std::ostream;

typedef std::complex<double> complex_t;

complex_t omega( int n )
{
    double arg = 2. * M_PI / n;
    return complex_t( cos(arg), sin(arg));
}


MatrixXcd omegaMatrix( int n )
{
    complex_t om = omega(n);
    MatrixXcd m;
    m.resize( n, n);

    int j,k;
    for( j = 0; j < n; j++)
        for( k = 0; k < n; k++)
            m(j,k) = pow( om, -k * j );
    
    return m / n;
}


VectorXcd DFT( VectorXcd Y )
{
    int n = Y.size();
    int start = (int)ceil( ((double)n) / 2. ) - 1;
    complex_t oms = pow( omega(n), (double)start);
    
    for( int k = 0; k < n; k++)
        Y[k] *= pow( oms, (double)k);
    
    VectorXcd res;
    res = omegaMatrix( n ) * Y;
    
    for( int k = 0; k < res.size(); k++)
        res[k] *= pow( -1., k + start);
    
    return res;
}


void print_input( const VectorXcd &Y, ostream &os)
{
     int n = Y.size();

     for( int k = 0; k < n; k++)
         os << (-M_PI + 2. * k * M_PI / n) << "  " << Y[k].real() << "\n";
}


void print_DFTPoly( const VectorXcd &Y, double step_x, ostream &os)
{
    int n = Y.size();
    
    VectorXcd coeffs = DFT( Y );

    int zero = (n - 1) / 2;
    
    for( double x = -M_PI; x <= M_PI; x += step_x)
    {
        double S = coeffs[zero].real();
        for( int k = 1; k < (int)ceil( ((double)n) / 2.); k++)
        {
            S += cos(k * x) * (coeffs[zero-k] + coeffs[zero+k]).real();
            S += sin(k * x) * (coeffs[zero-k] - coeffs[zero+k]).imag();
        }
        
        if( n % 2 == 0 )
            S += cos(n * x / 2.) * coeffs[coeffs.size()-1].real();
        
        os << x << "  " << S << "\n";
    }
}


int main( int argc, char **argv)
{
    VectorXcd in;

    int n = argc - 1;

    if( n < 3 )
    {
        std::cerr << "input too small\n";
        return 1;
    }
    
    in.resize(n);

    for( int i = 0; i < n; i++)
        in(i) = atof( argv[i+1] );
    
    ofstream plot_in( "plot_in.txt" );
    ofstream plot_out( "plot_out.txt" );
    print_input( in, plot_in);
    print_DFTPoly( in, 0.01, plot_out);
    
    return 0;
}
