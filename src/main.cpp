// By Cristiano Monteiro <cristianomonteiro@gmail.com> - 20220313
// Based on:
// http://numbers.computation.free.fr/Constants/Algorithms/nthdigit.html
// http://numbers.computation.free.fr/Constants/Algorithms/pidec.cpp

#include <Arduino.h>
#include <TM1637.h>
//#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

// ModInt should be a 64-bit long integer. It could be a double floating point type, provided
// adaptations of MulMod and SumMulMod functions are done
typedef int64_t ModInt;
ModInt _m;
double _invm;
long n;
String pd = "", last = "";

int CLK = 2;
int DIO = 3;

TM1637 tm(CLK, DIO);

int8_t PiDisp[] = {0x00, 0x00, 0x00, 0x00};

void InitializeModulo(ModInt m)
{
  _m = m;
  _invm = 1. / (double)m;
}

// Compute a*b modulo _m
inline ModInt MulMod(ModInt a, ModInt b)
{
  // classical trick to bypass the 64-bit limitation, when a*b does not fit into the ModInt type.
  // Works whenever a*b/_m is less than 2^52 (double type maximal precision)
  ModInt q = (ModInt)(_invm * (double)a * (double)b);
  return a * b - q * _m;
}

// Compute a*b+c*d modulo _m
inline ModInt SumMulMod(ModInt a, ModInt b, ModInt c, ModInt d)
{
  ModInt q = (ModInt)(_invm * ((double)a * (double)b + (double)c * (double)d));
  return a * b + c * d - q * _m;
}

// double MyTime()
// {
//   return ((double) clock())/ CLOCKS_PER_SEC;
// }

double FullDouble = 1024. * 1024. * 1024. * 1024. * 1024. * 8.; // 2^53

inline double easyround(double x)
{
  double y = x + FullDouble;
  y -= FullDouble;
  return y;
}

/* return g, A such that g=gcd(a,_m) and a*A=g mod _m  */
ModInt ExtendedGcd(ModInt a, ModInt &A)
{
  ModInt A0 = 1, A1 = 0;
  ModInt r0 = a, r1 = _m;

  while (r1 > 0.)
  {
    ModInt q = r0 / r1;

    ModInt tmp = A0 - q * A1;
    A0 = A1;
    A1 = tmp;

    tmp = r0 - q * r1;
    r0 = r1;
    r1 = tmp;
  }
  A = A0;
  return r0;
}

ModInt InvMod(ModInt a)
{
  ModInt A;
  a = a % _m;
  if (a < 0)
    a += _m;
  ModInt gcd = ExtendedGcd(a, A);
  // if (gcd!=1)
  // printf("pb, gcd should be 1\n");
  return A;
}

ModInt PowMod(ModInt a, long b)
{
  ModInt r, aa;

  r = 1;
  aa = a;
  while (1)
  {
    if (b & 1)
      r = MulMod(r, aa);
    b >>= 1;
    if (b == 0)
      break;
    aa = MulMod(aa, aa);
  }
  return r;
}

/* Compute sum_{j=0}^k binomial(n,j) mod m */
ModInt SumBinomialMod(long n, long k)
{
  // Optimisation : when k>n/2 we use the relation
  // sum_{j=0}^k binomial(n,j) =  2^n - sum_{j=0}^{n-k-1} binomial(n,j)
  //
  // Note : additionnal optimization, not afforded here, could be done when k is near n/2
  // using the identity sum_{j=0}^{n/2} = 2^(n-1) + 1/2 binomial(n,n/2). A global saving of
  // 20% or 25% could be obtained.
  if (k > n / 2)
  {
    ModInt s = PowMod(2, n) - SumBinomialMod(n, n - k - 1);
    if (s < 0)
      s += _m;
    return s;
  }
  //
  // Compute prime factors of _m which are smaller than k
  //
  const long NbMaxFactors = 20; // no more than 20 different prime factors for numbers <2^64
  long PrimeFactor[NbMaxFactors];
  long NbPrimeFactors = 0;
  ModInt mm = _m;
  // _m is odd, thus has only odd prime factors
  for (ModInt p = 3; p * p <= mm; p += 2)
  {
    if (mm % p == 0)
    {
      mm = mm / p;
      if (p <= k) // only prime factors <=k are needed
        PrimeFactor[NbPrimeFactors++] = p;
      while (mm % p == 0)
        mm = mm / p; // remove all powers of p in mm
    }
  }
  // last factor : if mm is not 1, mm is necessarily prime
  if (mm > 1 && mm <= k)
  {
    PrimeFactor[NbPrimeFactors++] = mm;
  }

  // BinomialExponent[i] will contain the power of PrimeFactor[i] in binomial
  long BinomialPower[NbMaxFactors];
  // NextDenom[i] and NextNum[i] will contain next multiples of PrimeFactor[i]
  long NextDenom[NbMaxFactors], NextNum[NbMaxFactors];
  for (long i = 0; i < NbPrimeFactors; i++)
  {
    BinomialPower[i] = 1;
    NextDenom[i] = PrimeFactor[i];
    NextNum[i] = PrimeFactor[i] * (n / PrimeFactor[i]);
  }

  ModInt BinomialNum0 = 1, BinomialDenom = 1;
  ModInt SumNum = 1;
  ModInt BinomialSecondary = 1;

  for (long j = 1; j <= k; j++)
  {
    // new binomial : b(n,j) = b(n,j-1) * (n-j+1) / j
    ModInt num = n - j + 1;
    ModInt denom = j;
    int BinomialSecondaryUpdate = 0;

    for (long i = 0; i < NbPrimeFactors; i++)
    {
      long p = PrimeFactor[i];
      // Test if p is a prime factor of num0
      if (NextNum[i] == n - j + 1)
      {
        BinomialSecondaryUpdate = 1;
        NextNum[i] -= p;
        BinomialPower[i] *= p;
        num /= p;
        while (num % p == 0)
        {
          BinomialPower[i] *= p;
          num /= p;
        }
      }
      // Test if p is a prime factor of denom0
      if (NextDenom[i] == j)
      {
        BinomialSecondaryUpdate = 1;
        NextDenom[i] += p;
        BinomialPower[i] /= p;
        denom /= p;
        while (denom % p == 0)
        {
          BinomialPower[i] /= p;
          denom /= p;
        }
      }
    }

    if (BinomialSecondaryUpdate)
    {
      BinomialSecondary = BinomialPower[0];
      for (long i = 1; i < NbPrimeFactors; i++)
        BinomialSecondary = MulMod(BinomialSecondary, BinomialPower[i]);
    }

    BinomialNum0 = MulMod(BinomialNum0, num);
    BinomialDenom = MulMod(BinomialDenom, denom);

    if (BinomialSecondary != 1)
    {
      SumNum = SumMulMod(SumNum, denom, BinomialNum0, BinomialSecondary);
    }
    else
    {
      SumNum = MulMod(SumNum, denom) + BinomialNum0;
    }
  }
  SumNum = MulMod(SumNum, InvMod(BinomialDenom));
  return SumNum;
}

/* return fractionnal part of 10^n*(a/b) */
double DigitsOfFraction(long n, ModInt a, ModInt b)
{
  InitializeModulo(b);
  ModInt pow = PowMod(10, n);
  ModInt c = MulMod(pow, a);
  return (double)c / (double)b;
}

/* return fractionnal part of 10^n*S, where S=4*sum_{k=0}^{m-1} (-1)^k/(2*k+1). m is even */
double DigitsOfSeries(long n, ModInt m)
{
  double x = 0.;
  for (ModInt k = 0; k < m; k += 2)
  {
    x += DigitsOfFraction(n, 4, 2 * k + 1) - DigitsOfFraction(n, 4, 2 * k + 3);
    x = x - easyround(x);
  }
  return x;
}

double DigitsOfPi(long n)
{
  double logn = log((double)n);
  long M = 2 * (long)(3. * n / logn / logn / logn);               // M is even
  long N = 1 + (long)((n + 15.) * log(10.) / (1. + log(2. * M))); // n >= N
  N += N % 2;                                                     // N should be even
  ModInt mmax = (ModInt)M * (ModInt)N + (ModInt)N;
  // printf("Parameters : M=%ld, N=%ld, M*N+M=%.0lf\n",M,N,(double)mmax);
  // double st = MyTime();
  double x = DigitsOfSeries(n, mmax);
  // printf("Series time : %.2lf\n",MyTime()-st);
  for (long k = 0.; k < N; k++)
  {
    ModInt m = (ModInt)2 * (ModInt)M * (ModInt)N + (ModInt)2 * (ModInt)k + 1;
    InitializeModulo(m);
    ModInt s = SumBinomialMod(N, k);
    s = MulMod(s, PowMod(5, N));
    s = MulMod(s, PowMod(10, n - N)); // n-N is always positive
    s = MulMod(s, 4);
    x += (2 * (k % 2) - 1) * (double)s / (double)m; // 2*(k%2)-1 = (-1)^(k-1)
    x = x - floor(x);
  }
  return x;
}

void setup()
{
  // Algorithm will only calculate from digit n (n=50, but I added 4 more to keep screen fluid) onwards, need to hard code the first ones
  String FirstDigits = "31415926535897932384626433832795028841971693993751058209";

  Serial.begin(9600);
  delay(200);
  Serial.println("Pidec...");

  // put your setup code here, to run once:
  tm.init();

  // set brightness; 0-7
  tm.set(6);

  while (FirstDigits.length() > 4)
  {
    Serial.println(FirstDigits);
    // Crude char to int conversion
    PiDisp[0] = FirstDigits.charAt(0) - '0';
    PiDisp[1] = FirstDigits.charAt(1) - '0';
    PiDisp[2] = FirstDigits.charAt(2) - '0';
    PiDisp[3] = FirstDigits.charAt(3) - '0';

    tm.display(PiDisp);

    FirstDigits = FirstDigits.substring(1);
    delay(500);
  }

  // Serial.println("(http://numbers.computation.free.fr/Constants/constants.html for more details)\n");
  // Serial.println("Enter n : ");
  // scanf("%ld",&n);
  n = 51;
  // if (n<50) {
  //   Serial.println("Error, n should be bigger than 50. Please retry\n");
  //   exit(1);
  // }
  // double st=MyTime();
}

void loop()
{
  //long int t1 = millis();
  double x = DigitsOfPi(n);
  double pow = 1.e9;
  double y = x * pow;
  // To be able to know exactly the digits of pi at position n, the
  // value (pow*x) should be not too close to an integer
  // Serial.println(y);
  while (pow > 10 && (y - floor(y) < 0.05 || y - floor(y) > 0.95))
  {
    pow /= 10.;
    y = x * pow;
  }
  // Serial.println("Digits of pi after n-th decimal digit : %.0lf\n",floor(y));

  if (last.charAt(1) == '0')
  {
    pd = last.substring(1);
    Serial.println("achei");
  }
  else
  {
    pd = String(y, 0);
  }

  last = pd;

  Serial.println(pd);

  // Crude char to int conversion
  PiDisp[0] = pd.charAt(0) - '0';
  PiDisp[1] = pd.charAt(1) - '0';
  PiDisp[2] = pd.charAt(2) - '0';
  PiDisp[3] = pd.charAt(3) - '0';

  tm.display(PiDisp);

  // printf("Total time: %.2lf\n",MyTime()-st);
  n++;
  //Serial.print("Time: "); Serial.print(millis()-t1); Serial.println(" milliseconds");
}
