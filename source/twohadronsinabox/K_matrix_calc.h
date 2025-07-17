#ifndef K_MATRIX_CALC_H
#define K_MATRIX_CALC_H

#include "K_matrix_info.h"
#include "fit_forms.h"
#include "xml_handler.h"
#include <map>
#include <memory>

// ***********************************************************************
// *                                                                     *
// *   The important classes "KtildeMatrixCalculator" and                *
// *   "KtildeInverseCalculator" are defined in this file.  Objects of   *
// *   these classes compute the numerical values of all of the matrix   *
// *   elements of the K-matrix or its inverse for a given block in      *
// *   the block diagonal basis.                                         *
// *                                                                     *
// *   The Ktilde and Ktildeinverse matrices computed here are           *
// *   defined by                                                        *
// *                                                                     *
// *    Kinv[aL',bL] = (qcm_a/mref)^(L'+1/2) Ktildeinv[aL',bL]           *
// *                              * (qcm_b/mref)^(L+1/2).                *
// *                                                                     *
// *   Note: this definition of Ktilde and Ktildeinv differs from that   *
// *   of Eq. (19) in Nucl. Phys. B924, 477 (2017).  The original        *
// *   definition caused an undesirable dependence on the lattice size   *
// *   of Ktildeinv.  We fix this problem here.                          *
// *                                                                     *
// *   Construction by XML content:                                      *
// *                                                                     *
// *   <KtildeMatrix> or <KtildeMatrixInverse>                           *
// *                                                                     *
// *     <Element>                                                       *
// *       <KElementInfo>...</KElementInfo>                              *
// *       <FitForm>...</FitForm>                                        *
// *     </Element>                                                      *
// *      ......                                                         *
// *                                                                     *
// *      <DecayChannels>                                                *
// *        <DecayChannelInfo>...</DecayChannelInfo>                     *
// *          ... other channels infos ...                               *
// *        (Order matters: first <DecayChannelInfo> tag is channel 0,   *
// *         second <DecayChannelInfo> is channel 1, and so on.  In the  *
// *         K-matrix, channels are referred to using the index 0, 1,...)*
// *      </DecayChannels>                                               *
// *                                                                     *
// *     <StartingValues>...</StartingValues>  (optional)                *
// *                                                                     *
// *   </KtildeMatrix> or </KtildeMatrixInverse>                         *
// *                                                                     *
// *   where                                                             *
// *                                                                     *
// *   <KElementInfo>                                                    *
// *      <JTimesTwo>6</JTimesTwo>                                       *
// *      <KIndex>...</KIndex>                                           *
// *      <KIndex>...</KIndex>    ( must have 1 or 2 <KIndex> tags)      *
// *   </KElementInfo>                                                   *
// *                                                                     *
// *   and                                                               *
// *                                                                     *
// *   <KIndex><L>3</L><Sx2>1<Sx2><Chan>0</Chan></KIndex>  OR            *
// *   <KIndex>L(3) 2S(1) chan(0)</KIndex> (order matters in this form)  *
// *                                                                     *
// *   For each decay channel, some necessary information is needed.     *
// *   The XML format for each decay channel information is              *
// *                                                                     *
// *    <DecayChannelInfo>                                               *
// *      <Particle1Name>pion</Particle1Name>                            *
// *      <Spin1TimesTwo>0</Spin1TimesTwo>                               *
// *      <Identical/> (if identical, do not include tags below)         *
// *      <Particle2Name>eta</Particle2Name>                             *
// *      <Spin2TimesTwo>2</Spin2TimesTwo>                               *
// *      <IntrinsicParities>same</IntrinsicParities> (or "opposite")    *
// *    </DecayChannelInfo>                                              *
// *                                                                     *
// *   The tag for each <IntrinsicParities> must match the tags in all   *
// *   other decay channels or an exception is thrown.                   *
// *                                                                     *
// *   The <FitForm> tag currently may include one or both of the        *
// *   following (if both, the full form is a SUM of the two subforms):  *
// *                                                                     *
// *   <Polynomial>           or      <Polynomial>                       *
// *     <Degree>3</Degree>             <Powers>0 1 2</Powers>           *
// *   </Polynomial>                  </Polynomial>                      *
// *                                                                     *
// *      -- form on left includes all power 0,1,...Degree               *
// *      -- form on right includes only those powers actually listed    *
// *                                                                     *
// *   <SumOfPoles>           or         <SumOfPoles>                    *
// *     <NumberOfPoles>2</NumberOfPoles>  <PoleIndices>0 2</PoleIndices>*
// *   </SumOfPoles>                     </SumOfPoles>                   *
// *                                                                     *
// *      -- form on left includes all pole indices 0,1,...Number-1      *
// *      -- form on right includes only those pole indices listed       *
// *                                                                     *
// *   Note: "KtildeMatrixCalculator" can use both <Polynomial> and      *
// *   <SumOfPoles> in <FitForm>,  BUT "KtildeInverseCalculator" can     *
// *   ONLY use <Polynomial>.                                            *
// *                                                                     *
// *   To assign initial values to the fit parameters, add the tag       *
// *   below:                                                            *
// *                                                                     *
// *     <StartingValues>                                                *
// *       <KFitParamInfo>                                               *
// *         <PolynomialTerm>                                            *
// *            <Power>3</Power>                                         *
// *            <KElementInfo>...</KElementInfo>                         *
// *         </PolynomialTerm>                                           *
// *         <StartingValue>0.4534</StartingValue>                       *
// *       </KFitParamInfo>                                              *
// *                                                                     *
// *       <KFitParamInfo>                                               *
// *         <PoleEnergy>                                                *
// *            <Index>3</Index>                                         *
// *            <JTimesTwo>2</JTimesTwo>                                 *
// *         </PoleEnergy>                                               *
// *         <StartingValue>2.2</StartingValue>                          *
// *       </KFitParamInfo>                                              *
// *                                                                     *
// *       <KFitParamInfo>                                               *
// *         <PoleCoupling>                                              *
// *            <Index>3</Index>                                         *
// *            <JTimesTwo>2</JTimesTwo>                                 *
// *            <KIndex>...</KIndex>                                     *
// *         </PoleCoupling>                                             *
// *         <StartingValue>3.3</StartingValue>                          *
// *       </KFitParamInfo>                                              *
// *                                                                     *
// *         ....                                                        *
// *     </StartingValues>                                               *
// *                                                                     *
// *   Alternatively, an object of one of these classes can be           *
// *   constructed via                                                   *
// *                                                                     *
// *   list<pair<KElementInfo,Polynomial> > pelems;                      *
// *   list<pair<KElementInfo,SumOfPoles> > selems;                      *
// *   list<pair<KElementInfo,SumOfPolesPlusPolynomial> > spelems;       *
// *   KtildeMatrixCalculator KC(pelems,selemn,spelems);                 *
// *                                                                     *
// ***********************************************************************

class KtildeMatrixCalculator {

  std::map<KElementInfo, FitForm*> m_fit;
  std::vector<KFitParamInfo> m_paraminfo;
  std::map<KFitParamInfo, uint> m_paramindices;
  std::vector<double> m_kappa_params;
  std::vector<DecayChannelInfo> m_decay_infos;

  // prevent copying, no default

  KtildeMatrixCalculator(const KtildeMatrixCalculator& inK);
  KtildeMatrixCalculator& operator=(const KtildeMatrixCalculator& inK);

public:
  // Default constructor for cloning
  KtildeMatrixCalculator();

  KtildeMatrixCalculator(XMLHandler& xmlin, bool require_initvals = false);

  KtildeMatrixCalculator(
      const std::list<std::pair<KElementInfo, Polynomial>>& pelems,
      const std::list<std::pair<KElementInfo, SumOfPoles>>& selems,
      const std::list<std::pair<KElementInfo, SumOfPolesPlusPolynomial>>&
          spelems,
      const std::vector<DecayChannelInfo>& chans);

  ~KtildeMatrixCalculator();

  // Deep copy/clone method to create an identical object with new pointers
  std::unique_ptr<KtildeMatrixCalculator> clone() const;

  uint getNumberOfParameters() const;

  void setParameterValues(std::vector<double> kappa_params);

  const std::vector<KFitParamInfo>& getFitParameterInfos() const {
    return m_paraminfo;
  }

  int getParameterIndex(
      const KFitParamInfo& kinfo) const; // returns -1 if not found

  double getParameterValue(const KFitParamInfo& kinfo) const;

  const std::vector<double>& getParameterValues() const {
    return m_kappa_params;
  }

  std::set<KElementInfo> getElementInfos() const;

  uint getNumberOfDecayChannels() const { return m_decay_infos.size(); }

  DecayChannelInfo getDecayChannelInfo(uint channel_index) const {
    return m_decay_infos.at(channel_index);
  }

  const std::vector<DecayChannelInfo>& getDecayChannelInfos() const {
    return m_decay_infos;
  }

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

  double calculate(uint Jtimestwo, uint Lp, uint Sptimestwo, uint chanp, uint L,
                   uint Stimestwo, uint chan, double Ecm_over_mref) const;

  bool isZero(uint Jtimestwo, uint Lp, uint Sptimestwo, uint chanp, uint L,
              uint Stimestwo, uint chan) const;

private:
  void initialize(std::list<XMLHandler>& xelems);

  void initialize_a_fitform(const KElementInfo& kinfo, XMLHandler& xmlin);

  void initialize_starting_values(XMLHandler& xmlin);

  friend class KtildeInverseCalculator;
};

// *************************************************************************

class KtildeInverseCalculator {

  std::map<KElementInfo, FitForm*> m_fit;
  std::vector<KFitParamInfo> m_paraminfo;
  std::map<KFitParamInfo, uint> m_paramindices;
  std::vector<double> m_kappa_params;
  std::vector<DecayChannelInfo> m_decay_infos;

  // disallow copying, no default

  KtildeInverseCalculator(const KtildeInverseCalculator& inK);
  KtildeInverseCalculator& operator=(const KtildeInverseCalculator& inK);

public:
  // Default constructor for cloning
  KtildeInverseCalculator();

  KtildeInverseCalculator(XMLHandler& xmlin, bool require_initvals = false);

  KtildeInverseCalculator(
      const std::list<std::pair<KElementInfo, Polynomial>>& polyelems,
      const std::vector<DecayChannelInfo>& chans);

  ~KtildeInverseCalculator();

  // Deep copy/clone method to create an identical object with new pointers
  std::unique_ptr<KtildeInverseCalculator> clone() const;

  uint getNumberOfParameters() const;

  void setParameterValues(std::vector<double> kappa_params);

  const std::vector<KFitParamInfo>& getFitParameterInfos() const {
    return m_paraminfo;
  }

  int getParameterIndex(
      const KFitParamInfo& kinfo) const; // returns -1 if not found

  double getParameterValue(const KFitParamInfo& kinfo) const;

  const std::vector<double>& getParameterValues() const {
    return m_kappa_params;
  }

  std::set<KElementInfo> getElementInfos() const;

  uint getNumberOfDecayChannels() const { return m_decay_infos.size(); }

  DecayChannelInfo getDecayChannelInfo(uint channel_index) const {
    return m_decay_infos.at(channel_index);
  }

  const std::vector<DecayChannelInfo>& getDecayChannelInfos() const {
    return m_decay_infos;
  }

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

  double calculate(uint Jtimestwo, uint Lp, uint Sptimestwo, uint chanp, uint L,
                   uint Stimestwo, uint chan, double Ecm_over_mref) const;

  bool isZero(uint Jtimestwo, uint Lp, uint Sptimestwo, uint chanp, uint L,
              uint Stimestwo, uint chan) const;

private:
  void initialize(std::list<XMLHandler>& xelems);

  void initialize_a_fitform(const KElementInfo& kinfo, XMLHandler& xmlin);

  void initialize_starting_values(XMLHandler& xmlin);
};

#endif
