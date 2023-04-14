
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {
  
  // draw raw value relative to a
  double ret = R::rnorm(mean, sd) - a;
  
  // reflect off boundries at 0 and (b-a)
  if (ret < 0 || ret > (b-a)) {
    
    // use multiple reflections to bring into range [-(b-a), 2(b-a)]
    if (ret < -(b - a)) {
      int n_double_intervals = floor(-ret / (b - a)) / 2;
      ret += 2 * (b - a) * (n_double_intervals + 1);
    } else if (ret > 2*(b - a)) {
      int n_double_intervals = floor(ret / (b - a) - 1) / 2;
      ret -= 2 * (b - a) * (n_double_intervals + 1);
    }
    
    // use one more reflection to bring into range [0, (b-a)]
    if (ret < 0) {
      ret = -ret;
    }
    if (ret > (b-a)) {
      ret = 2*(b-a) - ret;
    }
  }
  
  // no longer relative to a
  ret += a;
  
  // don't let ret equal exactly a or b
  if (ret == a) {
    ret += UNDERFLO_DOUBLE;
  } else if (ret == b) {
    ret -= UNDERFLO_DOUBLE;
  }
  
  return ret;
}
