#ifndef NI_REG_H
#define NI_REG_H
#include <vector>

// ******************************************************************
// *                                                                *
// *   This class is useful for regularizing the                    *
// *   Rummukainen-Gottlieb-Luescher (RGL) shifted zeta             *
// *   functions which have simple poles at the non-interacting     *
// *   energies.  Recall that                                       *
// *                                                                *
// *       s = shift vector                                         *
// *       gam = boost factor from lab to cm frame                  *
// *       usq = u^2 = (L qstar/(2 Pi))^2 for L^3 spatial lattice   *
// *                                                                *
// *   or in more detail,                                           *
// *                                                                *
// *     Total momentum P = (2 Pi / L) d   d = vector of integers   *
// *     Energy E evaluated in lab frame from MC simulations        *
// *     Mass of particle 1 is m1, mass of particle 2 is m2         *
// *                                                                *
// *        Ecm = sqrt(E^2 - P^2)                                   *
// *        gam = E / Ecm                                           *
// *        s = (1 + (m1^2-m2^2)/Ecm^2) * d                         *
// *        qcm^2 = Ecm^2/4 - (m1^2+m2^2)/2                         *
// *                   + (m1^2-m2^2)^2/(4 Ecm^2)                    *
// *        u^2 = L^2 qcm^2 / (2 Pi)^2                              *
// *                                                                *
// *   In an object of this class, one sets the "s" vector and the  *
// *   "gam" factor, then the regularizing function multiplies      *
// *   by a suitable tanh() factor for each non-interacting pole    *
// *   between a specified minimum and maximum value of u^2.        *
// *   The regularizing function is obtained by the                 *
// *   "getValue(usq)" member.                                      *
// *                                                                *
// *   The regularizing function is                                 *
// *                                                                *
// *         prod_k  tanh( (usq-uupoles[k])*alpha[k] )              *
// *                                                                *
// *   where the uupoles[k] are all of the non-interacting u^2      *
// *   values between the requested minimum and maximum u^2 values, *
// *   and the weights alpha[k] are chosen as follows:              *
// *                                                                *
// ******************************************************************

class NonIntRegulator {
  double m_gam, m_usq_min, m_usq_max;
  std::vector<double> m_svec;

  std::vector<double> m_non_int_uupoles;
  std::vector<double> m_tanh_weights;

  // Prevent copying.

  NonIntRegulator(const NonIntRegulator&);
  NonIntRegulator& operator=(const NonIntRegulator&);

public:
  NonIntRegulator();

  NonIntRegulator(const std::vector<double>& s, double gam, double usq_min,
                  double usq_max);

  ~NonIntRegulator() {}

  void reset(const std::vector<double>& s, double gam, double usq_min,
             double usq_max);

  void reset_if_diff(const std::vector<double>& s, double gam, double usq_min,
                     double usq_max);

  void clear();

  double getGamma() const { return m_gam; }

  double getUsqMin() const { return m_usq_min; }

  double getUsqMax() const { return m_usq_max; }

  std::vector<double> getSVector() const { return m_svec; }

  double getValue(double usq) const;

private:
  void setup();

  void setup_poles_noshift();

  void setup_poles_shifted();

  void setup_weights();

  double calcWeight(int k) const;
};

// **************************************************************
#endif
