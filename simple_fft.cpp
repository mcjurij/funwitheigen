// A simple fft. Done with hints from Prof. Weitz
//
// Compile:
// g++ -g -o simple_fft  simple_fft.cpp
//
// Execute:
// ./dft  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
//
// Plot using gnuplot:
// plot "plot_out.txt" using 1:2 with lines lw 3, "plot_in.txt" using 1:2:(0.03) with circles fill solid

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cassert>

using std::cout;
using std::vector;
using std::pow;
using std::ofstream;
using std::ostream;

typedef std::complex<double> complex_t;

static unsigned bit_reverse( unsigned int in, unsigned int used_bits)
{
    int r = 0;
    unsigned int i, bit_len = used_bits-1;
    
    for( i = 0; i <= bit_len; i++)
        r |= (in & 1 << i) ? (1 << (bit_len - i)) : 0;
    
    return r;
}


static void unscramble( vector<complex_t> &v )
{
    int n = v.size();
    int b = n - 1;  // sets all bits when n is a power of 2
    int used_bits;

    vector<complex_t> t( n );
    
    for( used_bits = 0; b>0; used_bits++)   // does a log2(n)
        b >>= 1;
    // cout << "used bits = " << used_bits << "\n";
    
    for( int i = 0; i < n; i++)
        t[ bit_reverse( i, used_bits) ] = v[i];
    
    v = t;
}


static complex_t omega( int n )
{
    double arg = 2. * M_PI / n;
    return complex_t( cos(arg), sin(arg));
}


static void butterfly( vector<complex_t> &a, size_t s, size_t e)  // s and e are indexes
{
    int d = e-s;
    assert( d > 0 );
    int n = d+1;
    complex_t om = 1. / omega( n );
    int h = n/2;      // half
    int hs = s + h;   // half start
    
    for( int i = 0; i < h; i++)
    {
        complex_t au = a.at( s + i );   // a from upper half
        complex_t al = a.at( hs + i );  // a from lower half
        
        a.at( s + i ) = au + al;
        a.at( hs + i ) = pow( om, (double)i) * (au - al);
    }

    if( h > 1 )
    {
        butterfly( a, s, hs - 1);
        butterfly( a, hs, e);
    }
}


vector<complex_t> FFT( const vector<complex_t> &in )  // in's length must be a power of 2
{
    vector<complex_t> Y = in;
    int n = Y.size();
    double start = ceil( ((double)n) / 2. ) - 1.;
    complex_t oms = pow( omega(n), start);
    
    for( int k = 0; k < n; k++)
        Y[k] *= pow( oms, (double)k);

    butterfly( Y, 0, n-1);
    unscramble( Y );

    for( int k = 0; k < n; k++)
        Y[k] *= pow( -1., k + start) / (double)n;
    
    return Y;
}


void print_input( const vector<complex_t> &Y, ostream &os)
{
     int n = Y.size();

     for( int k = 0; k < n; k++)
         os << (-M_PI + 2. * k * M_PI / n) << "  " << Y[k].real() << "\n";
}


void print_DFTPoly( const vector<complex_t> &Y, double step_x, ostream &os)
{
    int n = Y.size();
    
    vector<complex_t> coeffs = FFT( Y );
    
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
            S += cos(n * x / 2.) * coeffs[ coeffs.size()-1 ].real();
        
        os << x << "  " << S << "\n";
    }
}


int main( int argc, char **argv)
{
    vector<complex_t> in;
    int l = argc - 1;

    if( l < 3 )
    {
        std::cerr << "input too small\n";
        return 1;
    }
    
    int bm = 1;
    int bpos = 0;
    while( l-bm > 0 )
    {
        bm <<= 1;
        bpos++;
    }

    int n = 1<<bpos;

    in.resize( n );
    for( int i = 0; i < l; i++)
        in[i] = complex_t( atof( argv[i+1] ), 0.);
    
    if( l < n )
    {
        cout << "Filling with zeros up to index " << n << "\n";
        for( int i = l; i < n; i++)
            in[i] = complex_t( 0., 0.);
    }
    
    ofstream plot_in( "plot_in.txt" );
    ofstream plot_out( "plot_out.txt" );
    print_input( in, plot_in);
    print_DFTPoly( in, 0.01, plot_out);
    
    return 0;
}
