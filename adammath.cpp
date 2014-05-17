//file cotains source code for simple math functions not included in standard
//libraries or armadillo

//calculates the factorial of  a number
unsigned long long factorial(unsigned long long n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
// calculate nPp
unsigned long long permute(long long n,long long p)
{
    unsigned long long x=n;
    if(n>p)
    {
        for(long long iii=n-1;iii>n-p;iii--)
        {
            x*=iii;
        }
    }
    else if(n<p)
    {
        x=0;
    }
    return x;
}

// calculates nCp
unsigned long long choose(unsigned long long n, unsigned long long p)
{

    return (n < p) ? 0 : permute(n,p)/(factorial(n-p));
}
//evaluates as a kronecker delta
short kdelta(int n,int m)
{
    return (n==m)? 1: 0;
}
