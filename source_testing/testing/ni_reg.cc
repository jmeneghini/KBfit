#include "ni_reg.h"
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
using namespace std;

// ***********************************************************

NonIntRegulator::NonIntRegulator(const std::vector<double>& s, double gam,
                                 double usq_min, double usq_max)
    : m_gam(gam), m_usq_min(usq_min), m_usq_max(usq_max), m_svec(s) {
  setup();
}

NonIntRegulator::NonIntRegulator()
    : m_gam(1.0), m_usq_min(0.0), m_usq_max(0.0),
      m_svec(vector<double>(3, 0.0)) {
  setup();
}

void NonIntRegulator::reset(const std::vector<double>& s, double gam,
                            double usq_min, double usq_max) {
  m_gam = gam;
  m_usq_min = usq_min;
  m_usq_max = usq_max;
  m_svec = s;
  setup();
}

void NonIntRegulator::reset_if_diff(const std::vector<double>& s, double gam,
                                    double usq_min, double usq_max) {
  const double eps = 1e-10;
  if ((std::abs(m_gam - gam) > eps) || (std::abs(m_usq_min - usq_min) > eps) ||
      (std::abs(m_usq_max - usq_max) > eps) ||
      (std::abs(m_svec[0] - s[0]) > eps) ||
      (std::abs(m_svec[1] - s[1]) > eps) ||
      (std::abs(m_svec[2] - s[2]) > eps)) {
    reset(s, gam, usq_min, usq_max);
  }
}

void NonIntRegulator::clear() {
  m_gam = 1.0;
  m_usq_min = m_usq_max = 0.0;
  m_svec = vector<double>(3, 0.0);
  m_non_int_uupoles.clear();
  m_tanh_weights.clear();
}

void NonIntRegulator::setup() {
  if (m_usq_min > m_usq_max) {
    clear();
    throw(std::invalid_argument("Invalid set up of NonIntRegulator"));
  }
  double ss =
      m_svec[0] * m_svec[0] + m_svec[1] * m_svec[1] + m_svec[2] * m_svec[2];
  if (ss < 1e-15) {
    setup_poles_noshift();
  } else {
    setup_poles_shifted();
  }
  setup_weights();

  for (int k = 0; k < int(m_non_int_uupoles.size()); ++k) {
    cout << "pole " << k << "  at u^2 = " << m_non_int_uupoles[k]
         << " with weight " << m_tanh_weights[k] << endl;
  }
}

void NonIntRegulator::setup_poles_noshift() {
  cout << "set up poles no shift" << endl;
  int nmax = ceil(sqrt(m_usq_max));
  cout << "nmax = " << nmax << endl;
  set<int> values;
  for (int nx = 0; nx <= nmax; ++nx)
    for (int ny = nx; ny <= nmax; ++ny)
      for (int nz = ny; nz <= nmax; ++nz) {
        cout << "trying " << nx << " " << ny << " " << nz << " ";
        int nsq = nx * nx + ny * ny + nz * nz;
        cout << " nsq = " << nsq << endl;
        if ((double(nsq) >= m_usq_min) && (double(nsq) <= m_usq_max))
          values.insert(nx * nx + ny * ny + nz * nz);
      }
  m_non_int_uupoles.resize(values.size());
  int count = 0;
  for (set<int>::iterator it = values.begin(); it != values.end(); ++it) {
    m_non_int_uupoles[count] = double(*it);
    ++count;
  }
}

void NonIntRegulator::setup_poles_shifted() {
  /* m_non_int_uupoles.clear();
   m_tanh_weights.clear();

   double b1=-m_svec[0]/(2.0*m_gam);
   double b2=-m_svec[1]/(2.0*m_gam);
   double b3=-m_svec[2]/(2.0*m_gam);
   double ss=m_svec[0]*m_svec[0]+m_svec[1]*m_svec[1]+m_svec[2]*m_svec[2];
   double alpha=4.0*m_gam*(1.0-m_gam)/ss;



   usq=uu;*/
}

/*
   // does sub-sum for |n1| + |n2| + |n3| = nntotal
   // class T must have member K T(n1,n2,n3) --> return type K
   // class K must must have += defined, K=0 must assign to zero
   // and abs(K) must be defined.

template <typename T, typename K>
inline void do_nn_subsum(K& res, T& summand, int nntotal)
{
 for (int n3abs=0;n3abs<=nntotal;n3abs++){
    int n3step=(n3abs!=0)?2*n3abs:1;
    for (int n2abs=0;n2abs<=(nntotal-n3abs);n2abs++){
       int n2step=(n2abs!=0)?2*n2abs:1;
       int n1abs=nntotal-n3abs-n2abs;
       int n1step=(n1abs!=0)?2*n1abs:1;
    for (int n1=-n1abs;n1<=n1abs;n1+=n1step)
    for (int n2=-n2abs;n2<=n2abs;n2+=n2step)
    for (int n3=-n3abs;n3<=n3abs;n3+=n3step)
       res+=summand(n1,n2,n3);}}
}
*/

/*
double NonIntRegulator::pole_value(double b1, double b2, double b3, double
alpha, int n1, int n2, int n3) const
{
 double zb=1.0+alpha*(n1*b1+n2*b2+n3*b3);
 double z1=n1+zb*b1;
 double z2=n2+zb*b2;
 double z3=n3+zb*b3;
 return z1*z1+z2*z2+z3*z3;
}





inline cmplx ZetaSummand1::operator()(int n1, int n2, int n3) const
{
 double z1=n1;
 double z2=n2;
 double z3=n3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.value(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}



   //   Summation routines.


   // does sub-sum for |n1| + |n2| + |n3| = nntotal
   // class T must have member K T(n1,n2,n3) --> return type K
   // class K must must have += defined, K=0 must assign to zero
   // and abs(K) must be defined.

template <typename T, typename K>
inline void do_nn_subsum(K& res, T& summand, int nntotal)
{
 for (int n3abs=0;n3abs<=nntotal;n3abs++){
    int n3step=(n3abs!=0)?2*n3abs:1;
    for (int n2abs=0;n2abs<=(nntotal-n3abs);n2abs++){
       int n2step=(n2abs!=0)?2*n2abs:1;
       int n1abs=nntotal-n3abs-n2abs;
       int n1step=(n1abs!=0)?2*n1abs:1;
    for (int n1=-n1abs;n1<=n1abs;n1+=n1step)
    for (int n2=-n2abs;n2<=n2abs;n2+=n2step)
    for (int n3=-n3abs;n3<=n3abs;n3+=n3step)
       res+=summand(n1,n2,n3);}}
}

    //  Does the summation over n1, n2, n3.  Sum is done in order
    //  of increasing  |n1|+|n2|+|n3|=nn.  All terms up to "nnmin"
    //  are included, then after that, sum continues until convergence
    //  attained based on requested tolerances.

template <typename T, typename K>
inline void do_sum(K& res, T& summand, int nnmin, double abstol, double reltol,
                   bool oexclude=false)
{
 res=(oexclude)?0:summand(0,0,0);
 for (int nn=1;nn<=nnmin;nn++){
    do_nn_subsum<T,K>(res,summand,nn);}
 const int nnmax=1024;
 K add;
 for (int nn=nnmin+1;nn<=nnmax;nn++){
    add=0;
    do_nn_subsum<T,K>(add,summand,nn);
    res+=add;
    double addmag=abs(add);
    double resmag=abs(res);
    if ((addmag<abstol)||(addmag<reltol*resmag)) return;}
 throw(std::runtime_error("sum failed to converge"));
}

*/

void NonIntRegulator::setup_weights() {
  int n = m_non_int_uupoles.size();
  m_tanh_weights.resize(n);
  if (n == 0)
    return;
  for (int k = 0; k < n; ++k) {
    m_tanh_weights[k] = calcWeight(k);
  }
}

// used by "setup_weights"

double NonIntRegulator::calcWeight(int k) const {
  int n = m_non_int_uupoles.size();
  if (n == 1) {
    return 2.0;
  }
  double delta = 0.0;
  if ((k > 0) && (k < (n - 1))) {
    double p = m_non_int_uupoles[k];
    double pm = m_non_int_uupoles[k - 1];
    double pp = m_non_int_uupoles[k + 1];
    delta = std::min(0.5 * std::abs(pp - p), 0.5 * std::abs(p - pm));
  } else if (k == 0) {
    double p = m_non_int_uupoles[0];
    double pp = m_non_int_uupoles[1];
    delta = 0.5 * std::abs(pp - p);
  } else {
    double p = m_non_int_uupoles[n - 1];
    double pm = m_non_int_uupoles[n - 2];
    delta = 0.5 * std::abs(p - pm);
  }
  return 2.0 / delta;
}

double NonIntRegulator::getValue(double usq) const {
  double result = 1.0;
  if ((usq >= m_usq_min) && (usq < m_usq_max)) {
    int n = m_non_int_uupoles.size();
    for (int k = 0; k < n; ++k)
      result *= tanh((usq - m_non_int_uupoles[k]) * m_tanh_weights[k]);
  } else {
    throw(std::invalid_argument("u^2 out of range in NonIntRegulator"));
  }
  return result;
}

// *********************************************************************
/*
   //   Summation routines.


   // does sub-sum for |n1| + |n2| + |n3| = nntotal
   // class T must have member K T(n1,n2,n3) --> return type K
   // class K must must have += defined, K=0 must assign to zero
   // and abs(K) must be defined.

template <typename T, typename K>
inline void do_nn_subsum(K& res, T& summand, int nntotal)
{
 for (int n3abs=0;n3abs<=nntotal;n3abs++){
    int n3step=(n3abs!=0)?2*n3abs:1;
    for (int n2abs=0;n2abs<=(nntotal-n3abs);n2abs++){
       int n2step=(n2abs!=0)?2*n2abs:1;
       int n1abs=nntotal-n3abs-n2abs;
       int n1step=(n1abs!=0)?2*n1abs:1;
    for (int n1=-n1abs;n1<=n1abs;n1+=n1step)
    for (int n2=-n2abs;n2<=n2abs;n2+=n2step)
    for (int n3=-n3abs;n3<=n3abs;n3+=n3step)
       res+=summand(n1,n2,n3);}}
}

    //  Does the summation over n1, n2, n3.  Sum is done in order
    //  of increasing  |n1|+|n2|+|n3|=nn.  All terms up to "nnmin"
    //  are included, then after that, sum continues until convergence
    //  attained based on requested tolerances.

template <typename T, typename K>
inline void do_sum(K& res, T& summand, int nnmin, double abstol, double reltol,
                   bool oexclude=false)
{
 res=(oexclude)?0:summand(0,0,0);
 for (int nn=1;nn<=nnmin;nn++){
    do_nn_subsum<T,K>(res,summand,nn);}
 const int nnmax=1024;
 K add;
 for (int nn=nnmin+1;nn<=nnmax;nn++){
    add=0;
    do_nn_subsum<T,K>(add,summand,nn);
    res+=add;
    double addmag=abs(add);
    double resmag=abs(res);
    if ((addmag<abstol)||(addmag<reltol*resmag)) return;}
 throw(std::runtime_error("sum failed to converge"));
}

    //  Simple summation: no convergence check

template <typename T, typename K>
inline void do_simple_sum(K& res, T& summand, int nmax,
                          int offset1=0, int offset2=0, int offset3=0)
{
 res=0;
 for (int n1=offset1-nmax;n1<=nmax+offset1;n1++)
 for (int n2=offset2-nmax;n2<=nmax+offset2;n2++)
 for (int n3=offset3-nmax;n3<=nmax+offset3;n3++)
    res+=summand(n1,n2,n3);
}

template <typename T, typename K>
inline void do_simple_sum_oexclude(K& res, T& summand, int nmax,
                                   int offset1=0, int offset2=0, int offset3=0)
{
 res=0;
 for (int n1=offset1-nmax;n1<=nmax+offset1;n1++)
 for (int n2=offset2-nmax;n2<=nmax+offset2;n2++)
 for (int n3=offset3-nmax;n3<=nmax+offset3;n3++)
    if ((n1!=0)&&(n2!=0)&&(n3!=0)) res+=summand(n1,n2,n3);
}

// ***********************************************************

   // Gaussian quadrature vector

template <typename T>
class GQvec
{
   static const int m_ngauss=32;
   vector<T> m_store;

 public:

   GQvec() : m_store(m_ngauss) {}

   GQvec(const GQvec<T>& in) : m_store(in.m_store) {}

   GQvec(int c) : m_store(m_ngauss,T(c)) {}

   GQvec<T>& operator=(const GQvec<T>& in)
    {m_store=in.m_store; return *this;}

   GQvec<T>& operator=(int c)
    {for (int k=0;k<m_ngauss;k++) m_store[k]+=T(c); return *this;}

   GQvec<T>& operator+=(const GQvec<T>& add)
    {for (int k=0;k<m_ngauss;k++) m_store[k]+=add.m_store[k]; return *this;}

   int npoints() const {return m_ngauss;}

   double abs() const
    {double res=0.0;
     for (int k=0;k<m_ngauss;k++){
        double xx=std::abs(m_store[k]);
        if (xx>res) res=xx;}
     return res;}

   const T& operator[](int k) const {return m_store[k];}

   T& operator[](int k) {return m_store[k];}

};

template <typename T>
double abs(const GQvec<T>& in)
{
 return in.abs();
}




class ZetaSummand1
{
    double usq;
    YLMpoly ylm;

 public:

    ZetaSummand1(int lval, int mval, double usq);
    cmplx operator()(int n1, int n2, int n3) const;
};


ZetaSummand1::ZetaSummand1(int lval, int mval, double uu) : ylm(lval,mval)
{
 usq=uu;
}

inline cmplx ZetaSummand1::operator()(int n1, int n2, int n3) const
{
 double z1=n1;
 double z2=n2;
 double z3=n3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.value(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}


class ZetaSummand1s
{
    double alpha,b1,b2,b3,usq;
    YLMpoly ylm;

 public:

    ZetaSummand1s(int lval, int mval, const vector<double>& s, double ss,
                  double gam, double usq);
    cmplx operator()(int n1, int n2, int n3) const;
};


ZetaSummand1s::ZetaSummand1s(int lval, int mval, const vector<double>& s,
                             double ss, double gam, double uu) : ylm(lval,mval)
{
 b1=-s[0]/(2.0*gam);
 b2=-s[1]/(2.0*gam);
 b3=-s[2]/(2.0*gam);
 alpha=4.0*gam*(1.0-gam)/ss;
 usq=uu;
}

inline cmplx ZetaSummand1s::operator()(int n1, int n2, int n3) const
{
 double zb=1.0+alpha*(n1*b1+n2*b2+n3*b3);
 double z1=n1+zb*b1;
 double z2=n2+zb*b2;
 double z3=n3+zb*b3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.value(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}


  // *****

class ZetaSummand1_re
{
    double usq;
    YLMpoly ylm;

 public:

    ZetaSummand1_re(int lval, int mval, double usq);
    double operator()(int n1, int n2, int n3) const;
};


ZetaSummand1_re::ZetaSummand1_re(int lval, int mval, double uu) : ylm(lval,mval)
{
 usq=uu;
}

inline double ZetaSummand1_re::operator()(int n1, int n2, int n3) const
{
 double z1=n1;
 double z2=n2;
 double z3=n3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.realpart(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}


class ZetaSummand1s_re
{
    double alpha,b1,b2,b3,usq;
    YLMpoly ylm;

 public:

    ZetaSummand1s_re(int lval, int mval, const vector<double>& s, double ss,
                     double gam, double usq);
    double operator()(int n1, int n2, int n3) const;
};


ZetaSummand1s_re::ZetaSummand1s_re(int lval, int mval, const vector<double>& s,
                                   double ss, double gam, double uu) :
ylm(lval,mval)
{
 b1=-s[0]/(2.0*gam);
 b2=-s[1]/(2.0*gam);
 b3=-s[2]/(2.0*gam);
 alpha=4.0*gam*(1.0-gam)/ss;
 usq=uu;
}

inline double ZetaSummand1s_re::operator()(int n1, int n2, int n3) const
{
 double zb=1.0+alpha*(n1*b1+n2*b2+n3*b3);
 double z1=n1+zb*b1;
 double z2=n2+zb*b2;
 double z3=n3+zb*b3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.realpart(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}



  // *****

class ZetaSummand1_im
{
    double usq;
    YLMpoly ylm;

 public:

    ZetaSummand1_im(int lval, int mval, double usq);
    double operator()(int n1, int n2, int n3) const;
};


ZetaSummand1_im::ZetaSummand1_im(int lval, int mval, double uu) : ylm(lval,mval)
{
 usq=uu;
}

inline double ZetaSummand1_im::operator()(int n1, int n2, int n3) const
{
 double z1=n1;
 double z2=n2;
 double z3=n3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.imagpart(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}


class ZetaSummand1s_im
{
    double alpha,b1,b2,b3,usq;
    YLMpoly ylm;

 public:

    ZetaSummand1s_im(int lval, int mval, const vector<double>& s, double ss,
                     double gam, double usq);
    double operator()(int n1, int n2, int n3) const;
};


ZetaSummand1s_im::ZetaSummand1s_im(int lval, int mval, const vector<double>& s,
                                   double ss, double gam, double uu) :
ylm(lval,mval)
{
 b1=-s[0]/(2.0*gam);
 b2=-s[1]/(2.0*gam);
 b3=-s[2]/(2.0*gam);
 alpha=4.0*gam*(1.0-gam)/ss;
 usq=uu;
}

inline double ZetaSummand1s_im::operator()(int n1, int n2, int n3) const
{
 double zb=1.0+alpha*(n1*b1+n2*b2+n3*b3);
 double z1=n1+zb*b1;
 double z2=n2+zb*b2;
 double z3=n3+zb*b3;
 double zz=z1*z1+z2*z2+z3*z3;
 return ylm.imagpart(z1,z2,z3)/(zz-usq)*exp(-(zz-usq));
}


// ***************************************************************

class ZetaSummand2
{
    double usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<cmplx> summands;

 public:

    ZetaSummand2(int lval, int mval, double usq, const GQvec<double>& tvalues);
    const GQvec<cmplx>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2::ZetaSummand2(int lval, int mval, double uu,
                           const GQvec<double>& tvalues)
                     : ylm(lval,mval), tvals(tvalues)
{
 usq=uu;
}

inline const GQvec<cmplx>& ZetaSummand2::operator()(int n1, int n2, int n3)
const
{
 const double pipi=9.86960440108935861883447;
 double w1=n1;
 double w2=n2;
 double w3=n3;
 double ww=w1*w1+w2*w2+w3*w3;
 cmplx ylmval=ylm.value(w1,w2,w3);
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=ylmval*exp(-pipi*ww/tvals[k]);
 return summands;
}


class ZetaSummand2s
{
    double alpha,s1,s2,s3,usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<cmplx> summands;

 public:

    ZetaSummand2s(int lval, int mval, const vector<double>& s, double ss,
                  double gam, double usq, const GQvec<double>& tvalues);
    const GQvec<cmplx>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2s::ZetaSummand2s(int lval, int mval, const vector<double>& s,
                            double ss, double gam, double uu,
                            const GQvec<double>& tvalues)
                           : ylm(lval,mval), tvals(tvalues)
{
 s1=s[0]; s2=s[1]; s3=s[2];
 alpha=(gam-1.0)/ss;
 usq=uu;
}

inline const GQvec<cmplx>& ZetaSummand2s::operator()(int n1, int n2, int n3)
const
{
 const double pi=3.14159265358979323846264;
 const double pipi=9.86960440108935861883447;
 double wb=n1*s1+n2*s2+n3*s3;
 cmplx phase(cos(pi*wb),sin(pi*wb));
 wb*=alpha;
 double w1=n1+wb*s1;
 double w2=n2+wb*s2;
 double w3=n3+wb*s3;
 double ww=w1*w1+w2*w2+w3*w3;
 phase*=ylm.value(w1,w2,w3);
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=phase*exp(-pipi*ww/tvals[k]);
 return summands;
}


   // *********


class ZetaSummand2_re
{
    double usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<double> summands;

 public:

    ZetaSummand2_re(int lval, int mval, double usq, const GQvec<double>&
tvalues); const GQvec<double>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2_re::ZetaSummand2_re(int lval, int mval, double uu,
                                 const GQvec<double>& tvalues)
                     : ylm(lval,mval), tvals(tvalues)
{
 usq=uu;
}

inline const GQvec<double>& ZetaSummand2_re::operator()(int n1, int n2, int n3)
const
{
 const double pipi=9.86960440108935861883447;
 double w1=n1;
 double w2=n2;
 double w3=n3;
 double ww=w1*w1+w2*w2+w3*w3;
 double ylmval=ylm.realpart(w1,w2,w3);
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=ylmval*exp(-pipi*ww/tvals[k]);
 return summands;
}


class ZetaSummand2s_re
{
    double alpha,s1,s2,s3,usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<double> summands;

 public:

    ZetaSummand2s_re(int lval, int mval, const vector<double>& s, double ss,
                     double gam, double usq, const GQvec<double>& tvalues);
    const GQvec<double>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2s_re::ZetaSummand2s_re(int lval, int mval, const vector<double>& s,
                                   double ss, double gam, double uu,
                                   const GQvec<double>& tvalues)
                           : ylm(lval,mval), tvals(tvalues)
{
 s1=s[0]; s2=s[1]; s3=s[2];
 alpha=(gam-1.0)/ss;
 usq=uu;
}

inline const GQvec<double>& ZetaSummand2s_re::operator()(int n1, int n2, int n3)
const
{
 const double pi=3.14159265358979323846264;
 const double pipi=9.86960440108935861883447;
 double wb=n1*s1+n2*s2+n3*s3;
 double phase_re=cos(pi*wb);
 double phase_im=sin(pi*wb);
 wb*=alpha;
 double w1=n1+wb*s1;
 double w2=n2+wb*s2;
 double w3=n3+wb*s3;
 double ww=w1*w1+w2*w2+w3*w3;
 cmplx ylmval=ylm.value(w1,w2,w3);
 double phase=phase_re*ylmval.real()-phase_im*ylmval.imag();
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=phase*exp(-pipi*ww/tvals[k]);
 return summands;
}

   // *********


class ZetaSummand2_im
{
    double usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<double> summands;

 public:

    ZetaSummand2_im(int lval, int mval, double usq, const GQvec<double>&
tvalues); const GQvec<double>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2_im::ZetaSummand2_im(int lval, int mval, double uu,
                                 const GQvec<double>& tvalues)
                     : ylm(lval,mval), tvals(tvalues)
{
 usq=uu;
}

inline const GQvec<double>& ZetaSummand2_im::operator()(int n1, int n2, int n3)
const
{
 const double pipi=9.86960440108935861883447;
 double w1=n1;
 double w2=n2;
 double w3=n3;
 double ww=w1*w1+w2*w2+w3*w3;
 double ylmval=ylm.imagpart(w1,w2,w3);
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=ylmval*exp(-pipi*ww/tvals[k]);
 return summands;
}


class ZetaSummand2s_im
{
    double alpha,s1,s2,s3,usq;
    YLMpoly ylm;
    GQvec<double> tvals;
    mutable GQvec<double> summands;

 public:

    ZetaSummand2s_im(int lval, int mval, const vector<double>& s, double ss,
                     double gam, double usq, const GQvec<double>& tvalues);
    const GQvec<double>& operator()(int n1, int n2, int n3) const;
};


ZetaSummand2s_im::ZetaSummand2s_im(int lval, int mval, const vector<double>& s,
                                   double ss, double gam, double uu,
                                   const GQvec<double>& tvalues)
                           : ylm(lval,mval), tvals(tvalues)
{
 s1=s[0]; s2=s[1]; s3=s[2];
 alpha=(gam-1.0)/ss;
 usq=uu;
}

inline const GQvec<double>& ZetaSummand2s_im::operator()(int n1, int n2, int n3)
const
{
 const double pi=3.14159265358979323846264;
 const double pipi=9.86960440108935861883447;
 double wb=n1*s1+n2*s2+n3*s3;
 double phase_re=cos(pi*wb);
 double phase_im=sin(pi*wb);
 wb*=alpha;
 double w1=n1+wb*s1;
 double w2=n2+wb*s2;
 double w3=n3+wb*s3;
 double ww=w1*w1+w2*w2+w3*w3;
 cmplx ylmval=ylm.value(w1,w2,w3);
 double phase=phase_im*ylmval.real()+phase_re*ylmval.imag();
 for (int k=0;k<tvals.npoints();k++)
    summands[k]=phase*exp(-pipi*ww/tvals[k]);
 return summands;
}

// ***************************************************************


cmplx zeta_noshift(int l, int m, double usq)
{
 if (zetaRGL_RestFrame_IsZero(l,m)) return cmplx(0.0,0.0);
 const double pi=3.14159265358979323846264;

 bool oexclude=((fabs(usq)>1e-12)||(l==0))?false:true;
 cmplx res1;
 ZetaSummand1 summand1(l,m,usq);
 int nnmin=ceil(sqrt(std::abs(usq)));
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1,cmplx>(res1,summand1,nnmin,abstol,reltol,oexclude);

 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2 summand2(l,m,usq,tvals);
 GQvec<cmplx> intgrand;
 nnmin=6;
 do_sum<ZetaSummand2,GQvec<cmplx> >(intgrand,summand2,nnmin,abstol,reltol,true);

 cmplx res2(0.0,0.0);
 for (int k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 cmplx I(0.0,1.0),Il(1.0,0.0);
 for (int k=1;k<=l;k++) Il*=I;
 res2*=Il*pow(pi,1.5+double(l));

 if (l==0){
    res2+=pi*lzero_func(usq);}

 return res1+res2;
}


cmplx zeta_shifted(int l, int m, const vector<double>& s, double ss,
                   double gam, double usq)
{
 const double pi=3.14159265358979323846264;

 cmplx res1;
 ZetaSummand1s summand1(l,m,s,ss,gam,usq);
 int nnmin=ceil(std::abs(s[0]));
 int nn=ceil(std::abs(s[1])); if (nn>nnmin) nnmin=nn;
 nn=ceil(std::abs(s[2])); if (nn>nnmin) nnmin=nn;
 nn=ceil(sqrt(std::abs(usq))); if (nn>nnmin) nnmin=nn;
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1s,cmplx>(res1,summand1,nnmin,abstol,reltol);

 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2s summand2(l,m,s,ss,gam,usq,tvals);
 GQvec<cmplx> intgrand;
 nnmin=6;
 do_sum<ZetaSummand2s,GQvec<cmplx>
>(intgrand,summand2,nnmin,abstol,reltol,true);

 cmplx res2(0.0,0.0);
 for (int k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 cmplx I(0.0,1.0),Il(1.0,0.0);
 for (int k=1;k<=l;k++) Il*=I;
 res2*=Il*gam*pow(pi,1.5+double(l));

 if (l==0){
    res2+=gam*pi*lzero_func(usq);}

 return res1+res2;
}


 // ***********



double high_zeta_noshift_re(int l, int m, double usq, double abstol, double
reltol)
{
 const double pi=3.14159265358979323846264;
 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2_re summand2(l,m,usq,tvals);
 GQvec<double> intgrand;
 int nnmin=6;
 do_sum<ZetaSummand2_re,GQvec<double>
>(intgrand,summand2,nnmin,abstol,reltol,true); double res2=0.0; for (int
k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 res2*=pow(pi,1.5+double(l));
 return res2;
}

double high_zeta_noshift_im(int l, int m, double usq, double abstol, double
reltol)
{
 const double pi=3.14159265358979323846264;
 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2_im summand2(l,m,usq,tvals);
 GQvec<double> intgrand;
 int nnmin=6;
 do_sum<ZetaSummand2_im,GQvec<double>
>(intgrand,summand2,nnmin,abstol,reltol,true); double res2=0.0; for (int
k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 res2*=pow(pi,1.5+double(l));
 return res2;
}


double high_zeta_shifted_re(int l, int m, const vector<double>& s, double ss,
                            double gam, double usq, double abstol, double
reltol)
{
 const double pi=3.14159265358979323846264;
 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2s_re summand2(l,m,s,ss,gam,usq,tvals);
 GQvec<double> intgrand;
 int nnmin=6;
 do_sum<ZetaSummand2s_re,GQvec<double>
>(intgrand,summand2,nnmin,abstol,reltol,true); double res2=0.0; for (int
k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 res2*=gam*pow(pi,1.5+double(l));
 return res2;
}

double high_zeta_shifted_im(int l, int m, const vector<double>& s, double ss,
                            double gam, double usq, double abstol, double
reltol)
{
 const double pi=3.14159265358979323846264;
 GQvec<double> tvals,wts;
 assign_gq_abscissae_weights(tvals,wts);
 ZetaSummand2s_im summand2(l,m,s,ss,gam,usq,tvals);
 GQvec<double> intgrand;
 int nnmin=6;
 do_sum<ZetaSummand2s_im,GQvec<double>
>(intgrand,summand2,nnmin,abstol,reltol,true); double res2=0.0; for (int
k=0;k<tvals.npoints();k++){
    res2+=wts[k]*intgrand[k]*exp(tvals[k]*usq)/pow(tvals[k],1.5+double(l));}
 res2*=gam*pow(pi,1.5+double(l));
 return res2;
}






double zeta_noshift_re(int l, int m, double usq)
{
 if (zetaRGL_RestFrame_IsZero(l,m)) return 0.0;
 const double pi=3.14159265358979323846264;
 bool oexclude=((fabs(usq)>1e-12)||(l==0))?false:true;
 double res1;
 ZetaSummand1_re summand1(l,m,usq);
 int nnmin=ceil(sqrt(std::abs(usq)));
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1_re,double>(res1,summand1,nnmin,abstol,reltol,oexclude);

 double res2=0.0;
 int p=l%4;   // apply I^L
 if (p==0)
    res2=high_zeta_noshift_re(l,m,usq,abstol,reltol);
 else if (p==1)
    res2=-high_zeta_noshift_im(l,m,usq,abstol,reltol);
 else if (p==2)
    res2=-high_zeta_noshift_re(l,m,usq,abstol,reltol);
 else
    res2=high_zeta_noshift_im(l,m,usq,abstol,reltol);
 if (l==0){
    res2+=pi*lzero_func(usq);}
 return res1+res2;
}


double zeta_shifted_re(int l, int m, const vector<double>& s, double ss,
                      double gam, double usq)
{
 const double pi=3.14159265358979323846264;
 double res1;
 ZetaSummand1s_re summand1(l,m,s,ss,gam,usq);
 int nnmin=ceil(std::abs(s[0]));
 int nn=ceil(std::abs(s[1])); if (nn>nnmin) nnmin=nn;
 nn=ceil(std::abs(s[2])); if (nn>nnmin) nnmin=nn;
 nn=ceil(sqrt(std::abs(usq))); if (nn>nnmin) nnmin=nn;
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1s_re,double>(res1,summand1,nnmin,abstol,reltol);

 double res2=0.0;
 int p=l%4;   // apply I^L
 if (p==0)
    res2=high_zeta_shifted_re(l,m,s,ss,gam,usq,abstol,reltol);
 else if (p==1)
    res2=-high_zeta_shifted_im(l,m,s,ss,gam,usq,abstol,reltol);
 else if (p==2)
    res2=-high_zeta_shifted_re(l,m,s,ss,gam,usq,abstol,reltol);
 else
    res2=high_zeta_shifted_im(l,m,s,ss,gam,usq,abstol,reltol);
 if (l==0){
    res2+=gam*pi*lzero_func(usq);}
 return res1+res2;
}


 // *******


double zeta_noshift_im(int l, int m, double usq)
{
 if (zetaRGL_RestFrame_IsZero(l,m)) return 0.0;
 bool oexclude=((fabs(usq)>1e-12)||(l==0))?false:true;
 double res1;
 ZetaSummand1_im summand1(l,m,usq);
 int nnmin=ceil(sqrt(std::abs(usq)));
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1_im,double>(res1,summand1,nnmin,abstol,reltol,oexclude);

 double res2=0.0;
 int p=l%4;   // apply I^L
 if (p==0)
    res2=high_zeta_noshift_im(l,m,usq,abstol,reltol);
 else if (p==1)
    res2=high_zeta_noshift_re(l,m,usq,abstol,reltol);
 else if (p==2)
    res2=-high_zeta_noshift_im(l,m,usq,abstol,reltol);
 else
    res2=-high_zeta_noshift_re(l,m,usq,abstol,reltol);
 return res1+res2;
}


double zeta_shifted_im(int l, int m, const vector<double>& s, double ss,
                       double gam, double usq)
{
 double res1;
 ZetaSummand1s_im summand1(l,m,s,ss,gam,usq);
 int nnmin=ceil(std::abs(s[0]));
 int nn=ceil(std::abs(s[1])); if (nn>nnmin) nnmin=nn;
 nn=ceil(std::abs(s[2])); if (nn>nnmin) nnmin=nn;
 nn=ceil(sqrt(std::abs(usq))); if (nn>nnmin) nnmin=nn;
 nnmin+=6;
 double abstol=1e-18, reltol=1e-12;
 do_sum<ZetaSummand1s_im,double>(res1,summand1,nnmin,abstol,reltol);

 double res2=0.0;
 int p=l%4;   // apply I^L
 if (p==0)
    res2=high_zeta_shifted_im(l,m,s,ss,gam,usq,abstol,reltol);
 else if (p==1)
    res2=high_zeta_shifted_re(l,m,s,ss,gam,usq,abstol,reltol);
 else if (p==2)
    res2=-high_zeta_shifted_im(l,m,s,ss,gam,usq,abstol,reltol);
 else
    res2=-high_zeta_shifted_re(l,m,s,ss,gam,usq,abstol,reltol);
 return res1+res2;
}



 // *******************************************************************


cmplx zetaRGL(int l, int m, const vector<double>& s, double gam, double usq)
{
 double ss=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
 if (ss<1e-15){
    return zeta_noshift(l,m,usq);}
 return zeta_shifted(l,m,s,ss,gam,usq);
}

double zetaRGL_re(int l, int m, const vector<double>& s, double gam, double usq)
{
 double ss=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
 if (ss<1e-15){
    return zeta_noshift_re(l,m,usq);}
 return zeta_shifted_re(l,m,s,ss,gam,usq);
}

double zetaRGL_im(int l, int m, const vector<double>& s, double gam, double usq)
{
 double ss=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
 if (ss<1e-15){
    return zeta_noshift_im(l,m,usq);}
 return zeta_shifted_im(l,m,s,ss,gam,usq);
}



// ********************************************************************


WZetaRGLCalculator::WZetaRGLCalculator(const std::vector<double>& s,
                                       double gam, double usq)
   : m_gam(gam), m_usq(usq), m_svec(s)  {}


WZetaRGLCalculator::WZetaRGLCalculator()
   : m_gam(1.0), m_usq(0.0), m_svec(vector<double>(3,0.0))  {}


void WZetaRGLCalculator::reset(const std::vector<double>& s,
                               double gam, double usq)
{
 m_wr_lm_values.clear();
 m_wi_lm_values.clear();
 m_gam=gam;
 m_usq=usq;
 m_svec=s;
}


void WZetaRGLCalculator::reset_if_diff(const std::vector<double>& s,
                                       double gam, double usq)
{
 const double eps=1e-10;
 if ((std::abs(m_gam-gam)>eps)||(std::abs(m_usq-usq)>eps)||
     (std::abs(m_svec[0]-s[0])>eps)||(std::abs(m_svec[1]-s[1])>eps)||
     (std::abs(m_svec[2]-s[2])>eps)){
    m_wr_lm_values.clear();
    m_wi_lm_values.clear();
    m_gam=gam;
    m_usq=usq;
    m_svec=s;}
}


void WZetaRGLCalculator::clear()
{
 m_wr_lm_values.clear();
 m_wi_lm_values.clear();
}


  //  evaluates  x^n  given x^2 and integer n>=0

cmplx WZetaRGLCalculator::uint_power_from_sq(double xsq, uint n)
{
 uint nn=n;
 bool flag=nn&0x1u;
 double rr,ri;
 if (flag){
    if (xsq>=0.0){
       rr=sqrt(xsq); ri=0.0;}
    else{
       rr=0.0; ri=sqrt(-xsq);}}
 else{
    rr=1.0; ri=0.0;}
 double v=xsq;
 double p=1.0;
 nn>>=1;
 while (nn){
    if (nn&0x1u) p*=v;
    v*=v; nn>>=1;}
 return cmplx(rr*p,ri*p);
}



cmplx& WZetaRGLCalculator::eval_wr_lm(uint l, uint m)
{
 uint lmcode=encode_lm(l,m);
 std::map<uint,cmplx>::iterator it=m_wr_lm_values.find(lmcode);
 if (it!=m_wr_lm_values.end()){
    //cout << "wr["<<l<<","<<m<<"] got from map"<<endl;
    return it->second;}
 double res=zetaRGL_re(l,m,m_svec,m_gam,m_usq);
 res*=0.17958712212516656169/m_gam;   // multiply by 1/(gam*pi^(3/2))
 cmplx zres(res/uint_power_from_sq(m_usq,l+1));
 std::pair<std::map<uint,cmplx>::iterator,bool> rt
     =m_wr_lm_values.insert(make_pair(lmcode,zres));
 if (!(rt.second))
    throw(std::runtime_error("Could not insert into map"));
 //cout << "wr["<<l<<","<<m<<"] computed "<<rt.first->second<<endl;
 return rt.first->second;
}


cmplx& WZetaRGLCalculator::eval_wi_lm(uint l, uint m)
{
 uint lmcode=encode_lm(l,m);
 std::map<uint,cmplx>::iterator it=m_wi_lm_values.find(lmcode);
 if (it!=m_wi_lm_values.end()){
    //cout << "wi["<<l<<","<<m<<"] got from map"<<endl;
    return it->second;}
 double res=zetaRGL_im(l,m,m_svec,m_gam,m_usq);
 res*=0.17958712212516656169/m_gam;   // multiply by 1/(gam*pi^(3/2))
 cmplx zres(res/uint_power_from_sq(m_usq,l+1));
 std::pair<std::map<uint,cmplx>::iterator,bool> rt
     =m_wi_lm_values.insert(make_pair(lmcode,zres));
 if (!(rt.second))
    throw(std::runtime_error("Could not insert into map"));
 //cout << "wi["<<l<<","<<m<<"] computed"<<rt.first->second<<endl;
 return rt.first->second;
}

*/
// **********************************************************************
