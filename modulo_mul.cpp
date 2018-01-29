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

static int base = 10000;  // changing this means you have to find new omegas (which is very difficult)
static int exp  = 4;      // dito

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


vector<mint> FFT( const vector<mint> &in, bool inverse=false)
{
    vector<mint> Y = in;
    int n = Y.size();

    butterfly( Y, 0, n-1, 0, inverse);
    unscramble( Y );
    
    return Y;
}


vector<mint> IFFT( const vector<mint> &in )
{
    vector<mint> Y = FFT( in, true);

    mint d = find_inverse( Y.size(), mod);
    assert( d != 0 );
    for( size_t i = 0; i < Y.size(); i++)
        Y[i] = mul_mod( Y[i], d, mod);
    
    return Y;
}


struct rc_omegas {
    int  n;
    mint mod;
    vector<mint> omegas;
};


rc_omegas omegas[2] = {
    { 32,  356665537, { 12459978, 250152976, 180225010, 280067532, 356665536 }},
    { 32,  356094209, { 12460007, 40729184, 132302520, 312423919, 356094208 }}
};


bool pick_omegas( int n )
{
    int t;
    int butterfly_depth;
    bool found = true;
    
    int depth = 0;

    rc_omegas rc_omega = omegas[1];  // trial and error
    assert( n == rc_omega.n );
    mod = rc_omega.mod;

    for( butterfly_depth = 0; butterfly_depth < rc_omega.omegas.size(); butterfly_depth++)
    {
        omega_array[ butterfly_depth ] = rc_omega.omegas[ butterfly_depth ];
    }
    
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
        found = true;
    else
        found = false;

    if( found )
        cout << "mod = " << mod << ", omega = " << omega_array[0] << "\n";
    return found;
}


vector<mint>  string_to_num( const string &ns )
{
    vector<mint> nv;
    int e;
    mint v;
    size_t i;
    string s = ns;
    reverse( s.begin(), s.end());
    
    for( i = 0; i < (s.size() - exp); i += exp)
    {
        size_t k;
        e = 1;
        v = 0;
        if( i + exp >= s.size() )
            break;
        for( k = 0; k < exp; k++)
        {
            if( !isdigit( s[i+k] ) )
                std::cerr << "Not a digit '" << s[i+k] << "'\n";
            v += ( s[i+k] - '0' ) * e;
            e *= 10;
        }
        
        nv.push_back( v );
    }

    if( i < s.size() )
    {
        e = 1;
        v = 0;
        for( ; i < s.size(); i++)
        {
            v += ( s[i] - '0' ) * e;
            e *= 10;
        }
        
        nv.push_back( v );
    }
    
    return nv;
}

typedef vector<mint> num_t;

num_t  add( const num_t &a, const num_t &b)
{
    size_t m = std::max( a.size(), b.size());
    int carry = 0;
    num_t r;

    r.resize( m );

    for( size_t i = 0; i < m; i++)
    {
        int s = (i < a.size() ? a[i] : 0 ) + (i < b.size() ? b[i] : 0) + carry;
        carry = 0;
            
        r[ i ] = s % base;
        carry = (s - r[i]) / base;
    }

    if( carry > 0 )
        r.push_back( carry );

    return r;
}


num_t mul( const num_t &a, const num_t &b)
{
    int carry = 0;
    num_t r;
    size_t pos = 0;

    for( ; pos < b.size(); pos++)
    {
        num_t p;
        p.resize( a.size()+pos );
        for( size_t i = 0; i < pos; i++)
        {
            p[i] = 0;
        }
        carry = 0;
        for( size_t i = pos; i < a.size()+pos; i++)
        {
            
            int s = a[i-pos] * b[pos] + carry;
            
            p[ i ] = s % base;
            carry = (s - p[i]) / base;
            //cout << " carry = " << carry << "\n";
        }
        
        if( carry > 0 )
            p.push_back( carry );

        r = add( r, p);
    }
 
    num_t tmp;

    pos = r.size()-1;
    while( r[pos] == 0 )
        pos--;
    while( pos > 0 )
        tmp.push_back( r[pos--] );
     tmp.push_back( r[pos] );
    reverse( begin(tmp), end(tmp));
    
    return tmp;
}


bool compare( const num_t &a, const num_t &b)
{
    size_t pos = 0;

    if( a.size() != b.size() )
        return false;
    
    for( ; pos < b.size(); pos++)
    {
        if( a[pos] != b[pos] )
            return false;
    }

    return true;
}


void print_coeff( const vector<mint> &v, const num_t &comp)
{
    int i;
    
    /*for( i = 0; i < v.size(); i++)
    {
        printf( "%d\n", v[ v.size() - 1 - i]);
    }
    */
    
    num_t r, tmp;
    r.push_back( 0 );

    int pos = 0;
    for( i = 0; i < v.size(); i++, pos += exp)
    {
        num_t c;

        c.resize( i+1 );
        c[i] = v[i];
        r = add( r, c);
    }

    pos = r.size()-1;
    while( r[pos] == 0 )
        pos--;
    while( pos > 0 )
        tmp.push_back( r[pos--] );
     tmp.push_back( r[pos] );
    std::reverse( begin(tmp), end(tmp));
    r = tmp;
    
    for( i = 0; i < r.size(); i++)
    {
        mint v = r[ r.size() - 1 - i];

        printf( "%04d", v);
    }
    cout << "\n";

    if( compare( comp, r) )
        cout << "EQUAL!\n";
    else
        cout << "DIFFER!\n";
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


void multiply( const vector<mint> &pv, vector<mint> &qv)
{
    vector<mint> P, Q;
    
    P = FFT( pv );
    Q = FFT( qv );
    
    vector<mint> R;
    R.resize( 32 );
    int i;
    for( i = 0; i < R.size(); i++)
        R[i] = mul_mod( P[i], Q[i], mod);
    
    vector<mint> poly_coeff = IFFT( R );

    num_t comp = mul( pv, qv);
    
    print_coeff( poly_coeff, comp);
}


int main( int argc, char **argv)
{
    if( argc != 3 )
    {
        std::cerr << "need 2 numbers\n";
        return 1;
    }
    
    string ps = argv[1];
    string qs = argv[2];

    vector<mint> pv = string_to_num( ps );
    vector<mint> qv = string_to_num( qs );
    int l = (int)std::max( pv.size(), qv.size());
    
    int bm = 1;
    int bpos = 0;
    while( l-bm > 0 )
    {
        bm <<= 1;
        bpos++;
    }

    int n = 1<<bpos;
    if( n > 16 )
    {
        std::cerr << "input too large\n";
        return 1;
    }
    else
        n = 32;

    cout << "n = " << n << "\n";
    
    pv = fill_vector( pv, n);
    qv = fill_vector( qv, n);
    
    if( !pick_omegas( n ) )
    {
        std::cerr << "no residual class found\n";
        return 1;
    }
  
    multiply( pv, qv);
    
    return 0;
}
