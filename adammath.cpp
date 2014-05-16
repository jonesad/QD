//file cotains source code for simple math functions not included in standard
//libraries or armadillo

//calculates the factorial of  a number
unsigned long long factorial(unsigned long long n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// calculates nCp
unsigned long long choose(unsigned long long n, unsigned long long p)
{

    return (n < p) ? 0 : factorial(n)/(factorial(p)*factorial(n-p));
}
//evaluates as a kronecker delta
short kdelta(int n,int m)
{
    return (n==m)? 1: 0;
}
