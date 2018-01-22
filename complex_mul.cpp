// See https://www.youtube.com/watch?v=G4XiNDprjXA
//
// Compile:
// g++ -g -o complex_mul complex_mul.cpp
//
// Execute:
// ./complex_mul
// then compare poly_coeff with the values @26:56

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cassert>

using std::cout;
using std::vector;
using std::pow;


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


static void butterfly( vector<complex_t> &a, size_t s, size_t e, bool inverse)  // s and e are indexes
{
    int d = e-s;
    assert( d > 0 );
    int n = d+1;
    complex_t om = inverse ? omega( n ) : 1. / omega( n );
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
        butterfly( a, s, hs - 1, inverse);
        butterfly( a, hs, e, inverse);
    }
}


vector<complex_t> FFT( const vector<complex_t> &in, bool inverse=false)  // in's length must be a power of 2
{
    vector<complex_t> Y = in;
    int n = Y.size();

    butterfly( Y, 0, n-1, inverse);
    unscramble( Y );
    
    return Y;
}


vector<complex_t> IFFT( const vector<complex_t> &in )  // in's length must be a power of 2
{
    vector<complex_t> Y = FFT( in, true);

    for( size_t i = 0; i < Y.size(); i++)
        Y[i] /= Y.size();
    
    return Y;
}


int main( int argc, char **argv)
{
    vector<complex_t> P, p, Q, q;

    p.resize( 8 );
    q.resize( 8 );
    
    p[0] = complex_t( 1., 0.);
    p[1] = complex_t( 2., 0.);
    p[2] = complex_t(-3., 0.);
    p[3] = complex_t( 4., 0.);
    p[4] = complex_t( 0., 0.);
    p[5] = complex_t( 0., 0.);
    p[6] = complex_t( 0., 0.);
    p[7] = complex_t( 0., 0.);

    q[0] = complex_t( 3., 0.);
    q[1] = complex_t(-5., 0.);
    q[2] = complex_t( 2., 0.);
    q[3] = complex_t(-1., 0.);
    q[4] = complex_t( 0., 0.);
    q[5] = complex_t( 0., 0.);
    q[6] = complex_t( 0., 0.);
    q[7] = complex_t( 0., 0.);
    
    P = FFT( p );
    Q = FFT( q );

    cout << "P =\n";
    for( int i = 0; i < P.size(); i++)
        cout << P[i] << "\n";
    cout << "Q =\n";
    for( int i = 0; i < Q.size(); i++)
        cout << Q[i] << "\n";

    vector<complex_t> R;
    R.resize( 8 );
    int i,j;
    for( i = 0; i < R.size(); i++)
        R[i] = P[i] * Q[i];
    
    cout << "R =\n";
    for( i = 0; i < R.size(); i++)
        cout << R[i] << "\n";

    
    vector<complex_t> poly_coeff = IFFT( R );
    
    cout << "poly_coeff =\n";
    for( i = 0; i < poly_coeff.size(); i++)
        cout << poly_coeff[i].real() << "\n";

    return 0;
}
