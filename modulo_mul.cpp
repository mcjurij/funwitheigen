// See https://www.youtube.com/watch?v=ytkcYkzN1oI by Prof Weitz
// 
// Compile:
// g++ -g -o modulo_mul modulo_mul.cpp
//
// Execute:
// ./modulo_mul 5431 7842
// 5431 times 7842 equals 42589902

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <algorithm>    // std::min

using std::cout;
using std::vector;
using std::string;

typedef uint64_t mint;

// helpers from Wikipedia
#ifdef DONT_USE_LONG_DOUBLE
uint64_t mul_mod( uint64_t a, uint64_t b, uint64_t m)
{
   long double x;
   uint64_t c;
   int64_t r;
   if( a >= m )
       a %= m;
   if( b >= m )
       b %= m;
   x = a;
   c = x * b / m;
   r = (int64_t)(a * b - c * m) % (int64_t)m;
   return r < 0 ? r + m : r;
}

#else

uint64_t mul_mod( uint64_t a, uint64_t b, uint64_t m)
{
   uint64_t d = 0, mp2 = m >> 1;
   int i;
   if (a >= m)
       a %= m;
   if (b >= m)
       b %= m;
   for (i = 0; i < 64; ++i)
   {
       d = (d > mp2) ? (d << 1) - m : d << 1;
       if (a & 0x8000000000000000ULL)
           d += b;
       if (d > m) d -= m;
       a <<= 1;
   }
   return d;
}
#endif


uint64_t pow_mod( uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 1;
    while( b > 0 )
    {
        if(b % 2 == 1)
            r = mul_mod(r, a, m);
        b = b >> 1;
        a = mul_mod(a, a, m);
    }
    return r;
}


static unsigned bit_reverse( unsigned int in, unsigned int used_bits)
{
    int r = 0;
    unsigned int i, bit_len = used_bits-1;
    
    for( i = 0; i <= bit_len; i++)
        r |= (in & 1 << i) ? (1 << (bit_len - i)) : 0;
    
    return r;
}


static void unscramble( vector<mint> &v )
{
    int n = v.size();
    int b = n - 1;  // sets all bits when n is a power of 2
    int used_bits;

    vector<mint> t( n );
    
    for( used_bits = 0; b>0; used_bits++)   // does a log2(n)
        b >>= 1;
    
    for( int i = 0; i < n; i++)
        t[ bit_reverse( i, used_bits) ] = v[i];
    
    v = t;
}


static const mint mod = 233;

static mint omega( int n )
{
    int om = 0;
    
    switch( n )
    {
        case 2:
            om = 232;
            break;
        case 4:
            om = 144;
            break;
        case 8:
            om = 12;
            break;
    }
    
    return om;
}


static mint inv_omega( int n )
{
    int iom = 0;
    
    switch( n )
    {
        case 2:
            iom = 232;
            break;
        case 4:
            iom = 89;
            break;
        case 8:
            iom = 136;
            break;
    }
    
    return iom;
}


static void butterfly( vector<mint> &a, size_t s, size_t e, bool inverse)  // s and e are indexes
{
    int d = e-s;
    assert( d > 0 );
    int n = d+1;
    mint om = inverse ? inv_omega( n ) : omega( n );
    int h = n/2;      // half
    int hs = s + h;   // half start

    for( int i = 0; i < h; i++)
    {
        mint au = a.at( s + i );   // a from upper half
        mint al = a.at( hs + i );  // a from lower half
        
        a.at( s + i ) = ((au%mod) + (al%mod));
        a.at( hs + i ) = mul_mod( pow_mod( om, i, mod), (mod + ((au%mod) - (al%mod))) % mod, mod);
    }

    if( h > 1 )
    {
        butterfly( a, s, hs - 1, inverse);
        butterfly( a, hs, e, inverse);
    }
}


vector<mint> FFT( const vector<mint> &in, bool inverse=false)  // in's length must be 8
{
    vector<mint> Y = in;
    int n = Y.size();

    butterfly( Y, 0, n-1, inverse);
    unscramble( Y );
    
    return Y;
}



static int find_inverse( int a, int n) // from Wikipedia
{
    int t= 0, newt = 1;
    int r = n, newr = a;
    
    while( newr != 0 )
    {
        int quotient = r / newr;
        int ht = t, hr = r;
        t = newt;
        newt = ht - quotient * newt;
        r = newr;
        newr = hr - quotient * newr;
    }
    
    if( r > 1 )
        return 0;  // :-/ :-( :-{
    if( t < 0 )
        t += n;
    return t;
}

vector<mint> IFFT( const vector<mint> &in )  // in's length must be 8
{
    vector<mint> Y = FFT( in, true);

    mint d = find_inverse( Y.size(), mod);
    for( size_t i = 0; i < Y.size(); i++)
        Y[i] = mul_mod( Y[i], d, mod);
    
    return Y;
}


vector<mint> mk_vector( const string &s, int n)
{
    int i, mlen = std::min( (int)s.length(), n);
    vector<mint> v;
    v.resize( n );
    
    for( i = 0; i < mlen; i++)
    {
        char c = s[ s.length() - 1 - i ];
        if( isdigit( c ) )
            v[i] = c - '0';
        else
            break;
    }
    
    for( ; i < n; i++)
        v[i] = 0;
    
    return v;
}


int main( int argc, char **argv)
{
    vector<mint> P, p, Q, q;

    string ps = argv[1];
    string pq = argv[2];
    
    P = FFT( mk_vector( ps, 8) );
    Q = FFT( mk_vector( pq, 8) );
    
    vector<mint> R;
    R.resize( 8 );
    int i;
    for( i = 0; i < R.size(); i++)
        R[i] = mul_mod( P[i], Q[i], mod);
    
    vector<mint> poly_coeff = IFFT( R );
    
    mint e = 1;
    mint r = 0;
    for( i = 0; i < poly_coeff.size(); i++)
    {
        r += poly_coeff[i] * e;
        e *= 10;
    }
    cout << r << "\n";
    return 0;
}
