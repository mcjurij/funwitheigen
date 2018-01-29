

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <algorithm>    // std::min

using namespace std;

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


static mint power_by_n( mint b, int n )
{
    mint h, r = 0;
    switch( n )
    {
        case 2:
            r = b*b;
            break;
        case 4:
            r = b*b*b*b;
            break;
        case 8:
            r = b*b*b*b * b*b*b*b;
            break;
        case 16:
            h = power_by_n( b, 8);
            r = h * h;
            break;
        case 32:
            h = power_by_n( b, 16);
            r = h * h;
            break;
        case 64:
            h = power_by_n( b, 32);
            r = h * h;
            break;
    }

    return r;
}


static mint mod = 0;

static mint omega_array[10];         // n = 2^10 max
static mint inv_omega_array[10];

static vector<mint> primes;

static bool is_primitive( mint om, int nth, mint mod)
{
    size_t i;

    int test_primes[10] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29};

    for( i = 0; i < 10; i++)
    {
        mint p = test_primes[i];

        if( p >= nth )
            break;
        
        if( (nth%p) == 0 )
        {
            if( pow_mod( om, nth/p, mod) == 1 )
                return false;
        }
    }

    return true;
}


bool pick_omegas( int n, mint min_mod)
{
    int t;
    int butterfly_depth;
    bool found = false;
    mint bmod;
    size_t pidx;
    for( pidx = 0; pidx < primes.size() && !found; pidx++)
    {
        bmod = primes[pidx];
        if( bmod < min_mod )
            continue;
        
        //mod = bmod*n + 1;
        mod = bmod;
        for( t = mod/100; t < mod; t++)
        {
            int tn = n;
            int om=t;
            butterfly_depth = 0;
            
            while( !found && pow_mod( om, tn, mod) == 1 && (butterfly_depth == 0 ? is_primitive( om, tn, mod) : true) )
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

    if( found )
    {
        cout << "n = " << n << ", mod = " << mod << ", omega = ";
        for( int i = 0; i< butterfly_depth; i++)
            cout << omega_array[i] << ", ";
        cout << "\n";
    }
    return found;
}

bool pick_omegas2( int n, mint min_mod)
{
    int t;
    int butterfly_depth;
    bool found = false;
    mint bmod;
    size_t pidx, start_pidx;
    int cnt;
    omega_array[ 0 ] = min_mod;
    
    for( start_pidx = 0; start_pidx < primes.size(); start_pidx++)
    {
        bmod = primes[start_pidx];
        if( bmod > min_mod )
            break;
    }
    
    //mod = bmod*n + 1;
    for( cnt = 0; cnt < 100; cnt++)
    {
        found = false;
        min_mod = omega_array[ 0 ] + 27;
        
        for( t = min_mod; t < min_mod*10 && !found; t++)
        {
            
            for( pidx = start_pidx; pidx < primes.size() && !found; pidx++)
            {
                bmod = primes[pidx];
                mod = bmod;
                
                butterfly_depth = 0;
                int tn = n;
                int om=t;
                while( !found && pow_mod( om, tn, mod) == 1 && (butterfly_depth == 0 ? is_primitive( om, tn, mod) : true) )
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
        
        if( found )
        {
            cout << "    { " << n << ",  " << mod << ", { ";
            for( int i = 0; i< butterfly_depth; i++)
                cout << omega_array[i] << ", ";
            cout << "}},\n";
        }
    }
    return found;
}


int main( int argc, char **argv)
{
    if( argc != 3 )
    {
        cerr << "primes file needed. minimum mod needed\n";
        return 1;
    }
    
    ifstream primes_file( argv[1] );
    
    for( string line; getline( primes_file, line ); )
    {
        if( line.length() > 0 )
        {
            mint p = atoi( line.c_str() );
            primes.push_back( p );
        }
    }

    mint min_mod = atoi( argv[2] );
    
    for( int n = 8; n <= 64; n<<=1)
    {
        if( !pick_omegas2( n, min_mod) )
        {
            cout << "no residual class found for n = " << n << "\n";
            
        }
    }
    
    return 0;
}
