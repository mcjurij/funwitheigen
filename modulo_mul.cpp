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


static int find_inverse( int a, int m) // from Wikipedia
{
    int t= 0, newt = 1;
    int r = m, newr = a;
    
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
        t += m;
    return t;
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


static mint mod = 0;

static mint omega_array[7];         // n = 2^7 max
static mint inv_omega_array[7];

static mint omega( int depth )
{
    return omega_array[ depth ];
}


static mint inv_omega( int depth )
{
    return inv_omega_array[ depth ];
}


static void butterfly( vector<mint> &a, size_t s, size_t e, int depth, bool inverse)  // s and e are indexes
{
    int d = e-s;
    assert( d > 0 );
    int n = d+1;
    mint om = !inverse ? inv_omega( depth ) : omega( depth );
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
        butterfly( a, s, hs - 1, depth + 1, inverse);
        butterfly( a, hs, e, depth + 1, inverse);
    }
}


vector<mint> FFT( const vector<mint> &in, bool inverse=false)  // in's length must be 8
{
    vector<mint> Y = in;
    int n = Y.size();

    butterfly( Y, 0, n-1, 0, inverse);
    unscramble( Y );
    
    return Y;
}


vector<mint> IFFT( const vector<mint> &in )  // in's length must be 8
{
    vector<mint> Y = FFT( in, true);

    mint d = find_inverse( Y.size(), mod);
    assert( d != 0 );
    for( size_t i = 0; i < Y.size(); i++)
        Y[i] = mul_mod( Y[i], d, mod);
    
    return Y;
}


bool pick_omegas( int n )
{
    int t;
    int butterfly_depth;
    bool found = false;
    for( mod = 1000; mod < 100000 && !found; mod++)
    {
        for( t = 233; t < mod; t++)
        {
            int tn = n;
            int om=t;
            butterfly_depth = 0;
            
            while( pow_mod( om, tn, mod) == 1 && !found )
            {
                omega_array[ butterfly_depth++ ] = om;
               
                om = pow_mod( om, 2, mod);
               
                if( tn > 2 )
                    tn >>= 1;
                else
                    found = true;
            }
            
            if( found && omega_array[ butterfly_depth-1 ] == mod-1 )
                break;
            else
                found = false;
        }
        
        if( found )
        {
            int depth = 0;
            for( depth = 0; depth < butterfly_depth; depth++)
            {
                int inv = find_inverse( omega_array[ depth ], mod);
                if( inv != 0 && inv != 1 )
                {
                    inv_omega_array[ depth ] = inv;
                }
                else
                {
                    found = false;
                    break;
                }
            }
            
            if( found && depth == butterfly_depth )
                break;
            else
                found = false;
        }
    }
    
    return found;
}


vector<mint> mk_vector( const string &s )
{
    int i, mlen = (int)s.length();
    vector<mint> v;
    v.resize( mlen );
    
    for( i = 0; i < mlen; i++)
    {
        char c = s[ s.length() - 1 - i ];
        if( isdigit( c ) )
            v[i] = c - '0';
        else
            break;
    }
    
    return v;
}


vector<mint> fill_vector( const vector<mint> &v, int n)
{
    vector<mint> r = v;
    size_t i = v.size();
    r.resize( n );
    
    for( ; i < n; i++)
        r[i] = 0;              // probably not needed?

    return r;
}


int main( int argc, char **argv)
{
    vector<mint> P, Q;

    if( argc != 3 )
    {
        std::cerr << "need 2 numbers\n";
        return 1;
    }
    
    string ps = argv[1];
    string qs = argv[2];

    vector<mint> pv = mk_vector( ps );
    vector<mint> qv = mk_vector( qs );
    int l = (int)std::max( pv.size(), qv.size());
    
    int bm = 1;
    int bpos = 0;
    while( l-bm > 0 )
    {
        bm <<= 1;
        bpos++;
    }

    int n = 1<< (bpos + 1);
    if( n > (1<<7) )
    {
        std::cerr << "input too large\n";
        return 1;
    }
    
    pv = fill_vector( pv, n);
    qv = fill_vector( qv, n);
    
    if( !pick_omegas( n ) )
    {
        std::cerr << "no residual class found\n";
        return 1;
    }
    
    P = FFT( pv );
    Q = FFT( qv );
    
    vector<mint> R;
    R.resize( n );
    int i;
    for( i = 0; i < R.size(); i++)
        R[i] = mul_mod( P[i], Q[i], mod);
    
    vector<mint> poly_coeff = IFFT( R );
    
    mint e = 1;
    mint r = 0;                   // 64 bit is too small to hold sum
    for( i = 0; i < poly_coeff.size(); i++)
    {
        r += poly_coeff[i] * e;   // you need something bigger here ;-) 
        e *= 10;
    }
    cout << r << "\n";
    return 0;
}
