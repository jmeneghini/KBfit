#ifndef BOX_QUANT_H
#define BOX_QUANT_H

#include "K_matrix_calc.h"
#include "box_matrix.h"
#include "cmframe.h"
#include "matrix.h"
#include "xml_handler.h"
#include <complex>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

typedef unsigned int uint;
typedef std::complex<double> cmplx;

// ***********************************************************************
// *                                                                     *
// *   The important class "BoxQuantization" is defined in this file.    *
// *   This class computes                                               *
// *                                                                     *
// *         det(Ktildeinverse-B)   or   det(1-Ktilde*B)                 *
// *                                                                     *
// *   where B is the all-important box matrix, and Ktilde is related    *
// *   to the scattering K-matrix through the effective range            *
// *   expansion                                                         *
// *                                                                     *
// *    Kinv[aL',bL] = (q_a/mref)^(L'+1/2) Ktildeinv[aL',bL]             *
// *                              * (q_b/mref)^(L+1/2).                  *
// *                                                                     *
// *   Note: this definition of Ktilde and Ktildeinv differs from that   *
// *   used in Nucl. Phys. B924, 477 (2017).  The original definition    *
// *   caused a dependence on the lattice size of Ktildeinv.  We need    *
// *   to avoid this here.                                               *
// *                                                                     *
// *   A separate object of this class is needed for each **block** for  *
// *   a single ensemble of Monte Carlo configurations.  A block is      *
// *   specified by total momentum P, little group irrep Lambda, and     *
// *   details of each channel: for each decay channel, we need to know  *
// *   the two particle spins, whether or not they are identical, and    *
// *   the product of their intrinsic parities.                          *
// *                                                                     *
// *   An object of this class can be constructed from input XML of the  *
// *   form                                                              *
// *                                                                     *
// *    <BoxQuantization>                                                *
// *      <TotalMomentumRay>oa</TotalMomentumRay>                        *
// *          ["ar"=at rest (0,0,0), "oa"=on axis (0,0,1),               *
// *           "pd"=planar diagonal (0,1,1),                             *
// *           "cd"=cubic diagonal (1,1,1)]                              *
// *      <TotalMomentumIntSquared>2</TotalMomentumIntSquared>           *
// *          [if total momentum is (2*Pi/L)(1,1,1), then IntSquared=3]  *
// *      <LGIrrep>T1u</LGIrrep>                                         *
// *        <LmaxValues>4 4</LmaxValues>  (one for each channel)         *
// *    </BoxQuantization>                                               *
// *                                                                     *
// *   The block is specified by a momentum ray (a string).  The         *
// *   momentum ray can take the values "ar" at rest P=0, "oa" on-axis   *
// *   P=(0,0,n), "pd" planar-diagonal P=(0,n,n), or "cd" cubic-diagonal *
// *   P=(n,n,n). The total momentum P is given in terms of three        *
// *   integers dx,dy,dz,  where P  = (2*Pi/L) * (dx,dy,dz), and "L"     *
// *   is the spatial length of the lattice                              *
// *                                                                     *
// *   The constructor also requires that the end user supply either     *
// *   a pointer to a KtildeMatrixCalculator object or a                 *
// *   KtildeInverseCalculator object.  If a KtildeMatrixCalculator is   *
// *   given, then  det(1-B*Ktilde) is computed.  If a                   *
// *   KtildeInverseCalculator pointer is given, then  det(Ktildeinv-B)  *
// *   is computed.                                                      *
// *                                                                     *
// *                                                                     *
// *   Once constructed, the most important tasks are to reset the       *
// *   decay particle masses, the length of the lattice, and to          *
// *   compute the quantization determinant for a given lab frame        *
// *   energy.  All energies and lengths are expressed in terms of       *
// *   some reference scale "mref".                                      *
// *                                                                     *
// *   To reset the length of the lattice and the decay particle         *
// *   masses of channel "a", where a=0,1,..., use member functions      *
// *                                                                     *
// *       setRefMassL(mref_L);                                          *
// *       setMassesOverRef(channel_index,mass1_over_ref,                *
// *                        mass2_over_ref);                             *
// *                                                                     *
// *   Some of the quantities that can be computed are below:            *
// *                                                                     *
// *      //  computes Omega(mu,Ktildeinv-B) or                          *
// *      //  Omega(mu,1-B*K)   (depending on K or Kinv mode)            *
// *                                                                     *
// *     double getOmegaFromElab(mu, Elab_over_mref);                    *
// *     double getOmegaFromEcm(mu, Ecm_over_mref);                      *
// *     ComplexHermitianMatrix B;                                       *
// *     RealSymmetricMatrix Kv;                                         *
// *     double getOmega(mu,Kv,B);                                       *
// *                                                                     *
// *      //  computes [det(Ktildeinv-B)]^(1/Ndet) or                    *
// *      //  [det(1-B*K)]^(1/Ndet) where Ndet is positive odd integer   *
// *      //    (depending on K or Kinv mode)                            *
// *                                                                     *
// *     double getDeterminantRootFromElab(Elab_over_mref, Ndet);        *
// *     double getDeterminantRootFromEcm(Ecm_over_mref, Ndet);          *
// *                                                                     *
// *   For the residual determinant fitting method, extract the          *
// *   box matrix at the observed energy, then the Ktilde or Ktildeinv   *
// *   separately, then there is a determinant routine that uses these   *
// *   extracted matrices:                                               *
// *                                                                     *
// *     ComplexHermitianMatrix B;                                       *
// *     RealSymmetricMatrix Kv;                                         *
// *     BQ.getBoxMatrixFromElab(Elab_over_mref,B);                      *
// *     BQ.getKtildeOrInverseFromElab(Elab_over_mref,Kv);               *
// *     double detroot=BQ.getDeterminantRoot(Kv,B,Ndet);                *
// *                                                                     *
// *      //  computes the real eigenvalues of Ktildeinv-B  or           *
// *      //  the eigenvalues of (B - B*K*B), divided by |detB|^(1/N)    *
// *      //  for NxN matrices (depending on K or Kinv mode)             *
// *                                                                     *
// *     vector<double> getEigenvaluesFromElab(Elab_over_mref);          *
// *     vector<double> getEigenvaluesFromEcm(Ecm_over_mref);            *
// *                                                                     *
// *                                                                     *
// *                                                                     *
// ***********************************************************************

//  Specifies a basis state for the Kinverse-B matrix in the quantization
//  condition:
//        | a S J L n >    a = channel index, n = occurrence

//   m_store encodes the channel index (6 bits),
//   total S times two (6 bits), total J times two (7 bits)
//   L (7 bits), occurrence index (6 bits)
//  Each state maintains a pointer to a BoxMatrix object.

class BoxQuantBasisState {
  BoxMatrix* m_boxmat;
  uint m_store;

  BoxQuantBasisState() = delete; // no default value

public:
  BoxQuantBasisState(BoxMatrix* boxmat, uint channel_index, uint Stimestwo,
                     uint Jtimestwo, uint L, uint occurrence);

  BoxQuantBasisState(const BoxQuantBasisState& in)
      : m_boxmat(in.m_boxmat), m_store(in.m_store) {}

  BoxQuantBasisState& operator=(const BoxQuantBasisState& in) {
    m_boxmat = in.m_boxmat;
    m_store = in.m_store;
    return *this;
  }

  ~BoxQuantBasisState() {}

  BoxMatrix* getBoxMatrixPtr() const { return m_boxmat; }

  uint getChannelIndex() const;

  uint getStimestwo() const;

  uint getJtimestwo() const;

  uint getL() const;

  uint getOccurrence() const;

  bool operator==(const BoxQuantBasisState& rhs) const;

  bool operator!=(const BoxQuantBasisState& rhs) const;

  bool operator<(const BoxQuantBasisState& rhs) const;

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

private:
  void check_encode(uint channel_index, uint Stimestwo, uint Jtimestwo, uint L,
                    uint occurrence);
};

inline uint BoxQuantBasisState::getChannelIndex() const {
  return m_store >> 26;
}

inline uint BoxQuantBasisState::getStimestwo() const {
  return (m_store >> 20) & 0x3Fu;
}

inline uint BoxQuantBasisState::getJtimestwo() const {
  return (m_store >> 13) & 0x7Fu;
}

inline uint BoxQuantBasisState::getL() const { return (m_store >> 6) & 0x7Fu; }

inline uint BoxQuantBasisState::getOccurrence() const {
  return m_store & 0x3Fu;
}

inline bool
BoxQuantBasisState::operator==(const BoxQuantBasisState& rhs) const {
  return (m_boxmat == rhs.m_boxmat) && (m_store == rhs.m_store);
}

inline bool
BoxQuantBasisState::operator!=(const BoxQuantBasisState& rhs) const {
  return (m_boxmat != rhs.m_boxmat) || (m_store != rhs.m_store);
}

inline bool BoxQuantBasisState::operator<(const BoxQuantBasisState& rhs) const {
  return ((m_boxmat < rhs.m_boxmat) ||
          ((m_boxmat == rhs.m_boxmat) && (m_store < rhs.m_store)));
}

// ******************************************************
// *                                                    *
// *   Declaration of the class "BoxQuantization".      *
// *                                                    *
// ******************************************************

//  Implementation details:

//  "m_basis" contains information about the basis of states.
//  "m_wzetas" contains pointers to the created "WZetaRGLCalculator"
//    objects; there is one for each channel since for a given
//    Elab, there is a different qcmsq,svec,usq for each channel.
//  "m_boxes" contains pointers to the "BoxMatrix" objects created;
//    there is one for each channel and each total S in each channel;
//    the elements of this list are a pair containing the pointer
//    to a "BoxMatrix" object and the index of the associated channel;
//    this channel index is needed for resetting masses in the
//    "BoxMatrix" objects.

class BoxQuantization {

  std::string m_lgirrep, m_lgirrepB;
  std::string m_momray;
  uint m_dx, m_dy, m_dz;
  std::vector<uint> m_Lmaxes;

  double m_mref_L;
  std::vector<double> m_masses1, m_masses2;

  std::set<BoxQuantBasisState> m_basis;
  std::list<std::pair<BoxMatrix*, uint>> m_boxes;
  std::list<WZetaRGLCalculator*> m_wzetas;

  KtildeMatrixCalculator* m_Kmat;
  KtildeInverseCalculator* m_Kinv;


  // Prevent copying and no default.

  BoxQuantization(const BoxQuantization&);
  BoxQuantization& operator=(const BoxQuantization&);
  BoxQuantization();

  bool is_box_matrix_inverse_root = true;
  double* Elab_min;
  double* Elab_max;

public:
  BoxQuantization(XMLHandler& xmlin, KtildeMatrixCalculator* Kmatptr);

  BoxQuantization(XMLHandler& xmlin, KtildeInverseCalculator* Kinvptr);

  BoxQuantization(
      XMLHandler& xmlin, KtildeMatrixCalculator* Kmatptr,
      KtildeInverseCalculator* Kinvptr); // one of the pointers must be null

  BoxQuantization(const std::string& mom_ray, uint mom_int_sq,
                  const std::string& lgirrep,
                  const std::vector<DecayChannelInfo>& chan_infos,
                  const std::vector<uint> lmaxes,
                  KtildeMatrixCalculator* Kmatptr);

  BoxQuantization(const std::string& mom_ray, uint mom_int_sq,
                  const std::string& lgirrep,
                  const std::vector<DecayChannelInfo>& chan_infos,
                  const std::vector<uint> lmaxes,
                  KtildeInverseCalculator* Kinvptr);

  ~BoxQuantization();

  struct EigenvalueRegularizingInfo {
    double E_min;
    double E_max;
    double in_scalar;
    double out_scalar;
  };

  std::string getMomRay() const { return m_momray; }

  std::string getLittleGroupIrrep() const { return m_lgirrep; }

  std::string getLittleGroupBoxIrrep() const { return m_lgirrepB; }

  uint getTotalMomentumIntegerSquared() const {
    return m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;
  }

  uint getNumberOfDecayChannels() const {
    return (m_Kmat) ? m_Kmat->getNumberOfDecayChannels()
                    : m_Kinv->getNumberOfDecayChannels();
  }

  DecayChannelInfo getDecayChannelInfo(uint channel_index) const {
    return (m_Kmat) ? m_Kmat->getDecayChannelInfo(channel_index)
                    : m_Kinv->getDecayChannelInfo(channel_index);
  }

  const std::vector<DecayChannelInfo>& getDecayChannelInfos() const {
    return (m_Kmat) ? m_Kmat->getDecayChannelInfos()
                    : m_Kinv->getDecayChannelInfos();
  }

  uint getLmax(uint channel_index) const { return m_Lmaxes.at(channel_index); }

  void setRefMassL(double mref_L);

  void setMassesOverRef(uint channel_index, double mass1_over_ref,
                        double mass2_over_ref);

  double getRefMassL() const { return m_mref_L; }

  double getMass1OverRef(uint channel_index) const {
    return m_masses1.at(channel_index);
  }

  double getMass2OverRef(uint channel_index) const {
    return m_masses2.at(channel_index);
  }

  bool isKtildeInverseMode() const { return (m_Kinv != 0); }

  bool isKtildeMode() const { return (m_Kmat != 0); }

  uint getNumberOfKtildeParameters() const {
    if (m_Kmat != 0)
      return m_Kmat->getNumberOfParameters();
    else
      return m_Kinv->getNumberOfParameters();
  }

  bool isBoxMatrixInverseRootMode() const {
    return is_box_matrix_inverse_root;
  }

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

  void outputBasis(XMLHandler& xmlout) const; // XML output

  std::string outputBasis(int indent = 0) const; // XML output

  uint getBasisSize() const { return m_basis.size(); }

  std::string getKeyString(int row, int col)
      const; // returns "B[momray,Psqint,IrrepB,2S,chan][2J',L',n'][2J,L,n]"

  void outputKFitParams(XMLHandler& xmlout) const; // XML output

  std::string outputKFitParams(int indent = 0) const; // XML output

  const std::vector<KFitParamInfo>& getFitParameterInfos() const;

  int getParameterIndex(
      const KFitParamInfo& kinfo) const; // returns -1 if not found

  double getParameterValue(const KFitParamInfo& kinfo) const;

  std::set<KElementInfo> getElementInfos() const;

  std::list<double> getFreeTwoParticleEnergies(double min_Elab_over_mref,
                                         double max_Elab_over_mref) const;

  double getEcmOverMrefFromElab(double Elab_over_mref) const;

  void getQcmsqOverMrefsqFromElab(double Elab_over_mref,
                                  RVector& qcmsq_over_mrefsq) const;

  void getBoxMatrixFromElab(double Elab_over_mref, ComplexHermitianMatrix& B);

  void getBoxMatrixFromElab(double Elab_over_mref, CMatrix& B);

  void getBoxMatrixFromEcm(double Ecm_over_mref, ComplexHermitianMatrix& B);

  void getBoxMatrixFromEcm(double Ecm_over_mref, CMatrix& B);

  void getKtildeFromElab(double Elab_over_mref, RealSymmetricMatrix& Ktilde);

  void getKtildeFromElab(double Elab_over_mref, RMatrix& Ktilde);

  void getKtildeFromEcm(double Ecm_over_mref, RealSymmetricMatrix& Ktilde);

  void getKtildeFromEcm(double Ecm_over_mref, RMatrix& Ktilde);

  void getKtildeinvFromElab(double Elab_over_mref,
                            RealSymmetricMatrix& Ktildeinv);

  void getKtildeinvFromElab(double Elab_over_mref, RMatrix& Ktildeinv);

  void getKtildeinvFromEcm(double Ecm_over_mref,
                           RealSymmetricMatrix& Ktildeinv);

  void getKtildeinvFromEcm(double Ecm_over_mref, RMatrix& Ktildeinv);

  // gets Ktilde or Ktildeinv, depending on mode
  void getKtildeOrInverseFromElab(double Elab_over_mref,
                                  RealSymmetricMatrix& KtildeOrInverse);

  void getKtildeOrInverseFromElab(double Elab_over_mref,
                                  RMatrix& KtildeOrInverse);

  void getKtildeOrInverseFromEcm(double Ecm_over_mref,
                                 RealSymmetricMatrix& KtildeOrInverse);

  void getKtildeOrInverseFromEcm(double Ecm_over_mref,
                                 RMatrix& KtildeOrInverse);

  void outputKBMatricesFromElab(double Elab_over_mref, std::ostream& fout);

  void outputKBMatricesFromEcm(double Ecm_over_mref, std::ostream& fout);

  //  computes Omega(mu,Ktildeinv-B) or Omega(mu,1-B*K)
  //         (depending on K or Kinv mode)

  double getOmegaFromElab(double mu, double Elab_over_mref);

  double getOmegaFromEcm(double mu, double Ecm_over_mref);

  double getOmega(double mu, const RealSymmetricMatrix& KtildeOrInverse,
                  const ComplexHermitianMatrix& B);

  //  computes [det(Ktildeinv-B)]^(1/Ndet) or [det(1-B*K)]^(1/Ndet)
  //  where Ndet is positive odd integer (depending on K or Kinv mode)

  double getDeterminantRootFromElab(double Elab_over_mref, uint Ndet);

  double getDeterminantRootFromEcm(double Ecm_over_mref, uint Ndet);

  double getDeterminantRoot(const RealSymmetricMatrix& KtildeOrInverse,
                            const ComplexHermitianMatrix& B, uint Ndet);

  //  computes the real eigenvalues of Ktildeinv-B  or
  //  the eigenvalues of (B - B*K*B), divided by |detB|^(1/N)
  //  for NxN matrices (depending on K or Kinv mode)

  std::vector<double> getEigenvaluesFromElab(double Elab_over_mref, EigenvalueRegularizingInfo* ev_reg_info = nullptr);

  std::vector<double> getEigenvaluesFromEcm(double Ecm_over_mref, EigenvalueRegularizingInfo* ev_reg_info = nullptr);

  //  computes [det(B)]^(1/Ndet)

  double getBoxMatrixDeterminantRootFromElab(double Elab_over_mref, uint Ndet);

  double getBoxMatrixDeterminantRootFromEcm(double Ecm_over_mref, uint Ndet);

  double getBoxMatrixDeterminantRoot(const ComplexHermitianMatrix& B,
                                     uint Ndet);

private:
  void clear();

  void set_dvector(uint mom_int_sq);

  void xmlinitialize(XMLHandler& xlmin);

  void initialize();

  void setup_basis();

  std::string getLGIrrepParityFlip();

  void assign(const cmplx& value, uint row, uint col, bool herm,
              ComplexHermitianMatrix& Bh, CMatrix& B);

  void assign(double value, uint row, uint col, bool herm,
              RealSymmetricMatrix& Kh, RMatrix& K);

  void assign_matrices(double E_over_mref, bool Elab, ComplexHermitianMatrix& B,
                       RealSymmetricMatrix& Kv);

  void get_box_matrix(double E_over_mref, ComplexHermitianMatrix& Bh,
                      CMatrix& B, bool Elab, bool herm);

  template <typename T>
  void get_ktilde_matrix(double E_over_mref, RealSymmetricMatrix& Kh,
                         RMatrix& K, bool Elab, bool herm, T* evalptr);

  void get_Vdiag_and_U(const ComplexHermitianMatrix& B, std::vector<double>& V_eigvals, CMatrix& U);

  template <typename T>
  std::set<BoxQuantBasisState> find_excluded_states_from_ktilde(T* evalptr);

  void output_matrices(double E_over_mref, bool Elab, std::ostream& fout);

  double get_determinant(double E_over_mref, bool Elab, uint Ndet);

  double get_determinant(uint N, const RealSymmetricMatrix& Kv,
                         const ComplexHermitianMatrix& B, uint Ndet);

  std::vector<double> get_eigenvalues(double E_over_mref, bool Elab, EigenvalueRegularizingInfo* ev_reg_info);

  double get_omega(double mu, double E_over_mref, bool Elab);

  double get_omega(double mu, uint N, const RealSymmetricMatrix& Kv,
                   const ComplexHermitianMatrix& B);
};

#endif
