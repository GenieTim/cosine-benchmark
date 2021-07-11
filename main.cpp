/**
 * Contributing Author: Tim Bernhard, ETHZ
 */

#include "timing_tsc_x86.h"
#include <cmath>
#include <iostream>
#include <random>

using namespace std;

#define NR_EXP 10000000

#define MY_PI 3.14159265358979323846  // pi
#define MY_2PI 6.28318530717958647692 // 2pi
#define MY_PI2 1.57079632679489661923 // pi/2
#define TWO_PI_INVERSE 0.1591549430918953357690176

/**
 * Approximation of cosine. |e(x)|≤2e-9
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical
 * functions: with formulas, graphs, and mathematical tables, p. 43
 */
static double abramowitzCosinePolynomial(double x) {
  const double x2 = x * x;

  const double term1 = x2 * -0.0000002605 + 0.0000247609;
  const double term2 = x2 * term1 - 0.0013888397;
  const double term3 = x2 * term2 + 0.0416666418;
  const double term4 = x2 * term3 - 0.4999999963;
  return 1.0 + x2 * term4;
}

/**
 * Approximation of cosine
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical
 * functions: with formulas, graphs, and mathematical tables, p. 43
 */
static double cosAbramowitz(double x) {
  // wrap x within [0, TWO_PI)
  const double a = x * TWO_PI_INVERSE;
  x -= static_cast<int>(a) * MY_2PI;
  if (x < 0.0f)
    x += MY_2PI;

  // 4 pieces of hills: wrap x within [0, pi/2]
  if (x < MY_PI2)
    return abramowitzCosinePolynomial(x);
  else if (x < MY_PI)
    return -abramowitzCosinePolynomial(MY_PI - x);
  else if (x < 3.0f * MY_PI2)
    return -abramowitzCosinePolynomial(x - MY_PI);
  else
    return abramowitzCosinePolynomial(MY_2PI - x);
}

static double sinKohlmeyer(double x) {
  int k;
  double y;
  double z;
  union {
    int i;
    double d;
  } di;

  z = x;
  z *= 0.3183098861837907;
  di.d = z + 6755399441055744.0;
  z = k = di.i;
  z *= 3.1415926535897932;
  x -= z;
  y = x;
  y *= x;
  z = 0.0073524681968701;
  z *= y;
  z -= 0.1652891139701474;
  z *= y;
  z += 0.9996919862959676;
  x *= z;
  k &= 1;
  k += k;
  z = k;
  z *= x;
  x -= z;

  return x;
}

static double cosKohlmeyer(double x) { return sinKohlmeyer(x + 0.5 * M_PI); }

static double cosStd(double x) { return cos(x); }

/**
 * Main entrypoint of benchmark
 */
int main(int argc, char const *argv[]) {

  auto angles = new double[NR_EXP];
  auto results = new double[NR_EXP];

  std::uniform_real_distribution<double> unif(0, 2 * MY_2PI);
  std::default_random_engine re;
  auto engine = re;
  engine.seed(5555);

  for (int i = 0; i < NR_EXP; ++i) {
    angles[i] = unif(engine);
    results[i] = cos(angles[i]);
  }

  myInt64 startTime;
  myInt64 endTime;
  double absoluteError;

#define RUN_STUDY(cosineFunction)                                              \
  absoluteError = 0;                                                           \
  startTime = start_tsc();                                                     \
  for (int i = 0; i < NR_EXP; i += 1) {                                        \
    absoluteError += (results[i] - cosineFunction(angles[i])) / NR_EXP;        \
  }                                                                            \
  endTime = stop_tsc(startTime);                                               \
  cout << endTime << ": " << absoluteError << ": " << #cosineFunction << endl;

  RUN_STUDY(cosStd)
  RUN_STUDY(cosKohlmeyer)
  RUN_STUDY(cosAbramowitz)

  delete[] angles;
  return 0;
}
