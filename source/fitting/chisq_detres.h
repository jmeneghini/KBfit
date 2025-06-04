#ifndef CHISQ_DETRES_H
#define CHISQ_DETRES_H

#include "box_quant.h"
#include "chisq_base.h"
#include "kbobs_handler.h"
#include "matrix.h"
#include "xml_handler.h"

// **************************************************************************
// *                                                                        *
// *    The class "DeterminantResidualFit", derived from the                *
// *    base class "ChiSquare", is defined here.  It evaluates the          *
// *    chi^2 value associated with fits using the Determinant              *
// *    Residual method (see Nucl. Phys. B924, 477 (2017)).  In the         *
// *    residuals, either of the following quantities is used:              *
// *                                                                        *
// *         Omega(mu, Ktildeinverse-B)   or   Omega(mu, 1-Ktilde*B)        *
// *                                                                        *
// *    where                                                               *
// *                                det(A)                                  *
// *         Omega(mu, A) = -----------------------                         *
// *                        det( (mu+A*A^T)^(1/2) )                         *
// *                                                                        *
// *    B is the all-important box matrix, and Ktilde is related to the     *
// *    scattering K-matrix through the effective range expansion           *
// *                                                                        *
// *     Kinv[aL',bL] = (q_a/mref)^(L'+1/2) Ktildeinv[aL',bL](Ecm/mref)     *
// *                               * (q_b/mref)^(L+1/2).                    *
// *                                                                        *
// *    Note: this definition of Ktilde and Ktildeinv differs from that     *
// *    used in Nucl. Phys. B924, 477 (2017).  The original definition      *
// *    caused a dependence on the lattice size of Ktildeinv.  We need      *
// *    to avoid this here.  The input XML for constructing such an         *
// *    object is shown below.                                              *
// *                                                                        *
// *                                                                        *
// *    <DeterminantResidualFit>                                            *
// *                                                                        *
// *      <OmegaMu>8.0</OmegaMu>  (optional)                                *
// *                                                                        *
// *      <Verbose/>  (optional)                                            *
// *                                                                        *
// *      <KtildeMatrixInverse>  (or <KtildeMatrix>)                        *
// *          .....                                                         *
// *      </KtildeMatrixInverse>                                            *
// *                                                                        *
// *      <DefaultEnergyFormat>reference_ratio</DefaultEnergyFormat>        *
// *                   (or time_spacing_product)                            *
// *                                                                        *
// *      <MCEnsembleParameters>...</MCEnsembleParameters>                  *
// *        ... one for each Monte Carlo ensemble                           *
// *                                                                        *
// *      <KBBlock>...</KBBlock>                                            *
// *        ... one for each KB quantization block                          *
// *                                                                        *
// *      <KBObservables>                                                   *
// *       <MCSamplingInfo>...</MCSamplingInfo>                             *
// *       <SamplingData>                                                   *
// *          <FileName>...</FileName>  (all sampling files needed to       *
// *             ....                    obtain all above MCObsInfo's)      *
// *       </SamplingData>                                                  *
// *      </KBObservables>                                                  *
// *                                                                        *
// *    </DeterminantResidualFit>                                           *
// *                                                                        *
// *                                                                        *
// *    Notes:                                                              *
// *                                                                        *
// *    -- An L^3 spatial lattice is required.  Eventually, the length      *
// *       times the reference scale is needed.  To determine this,         *
// *       the reference scale times the temporal lattice spacing must      *
// *       be specified in <ReferenceMassTimeSpacingProduct>.               *
// *                                                                        *
// *    -- If <OmegaMu> is specified, then the Omega function is used       *
// *       in the residual.  If absent, the determinant itself is used.     *
// *                                                                        *
// *    -- Either <KtildeMatrixInverse> or <KtildeMatrix> must be given.    *
// *       Depending on which is input, det(1-K*B) or det(K^(-1)-B)         *
// *       is used.                                                         *
// *                                                                        *
// *    -- If using an anisotropic lattice, a tag <LatticeAnisotropy>,      *
// *       which is the spatial over the temporal spacing, must be given.   *
// *       If this tag is absent, an isotropic lattice is assumed.          *
// *                                                                        *
// *    -- All lab-frame energies and particle masses can be input either   *
// *       as ratios of a reference mass or as a product with the time      *
// *       spacing of the lattice.  The default format should be specified  *
// *       in the tag <DefaultEnergyFormat> whose value can be either       *
// *       "reference_ratio" or "time_spacing_product".                     *
// *                                                                        *
// *                                                                        *
// *    For each Monte Carlo ensemble involved, there should be a tag       *
// *    with the format:                                                    *
// *                                                                        *
// *      <MCEnsembleParameters>                                            *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <ReferenceMassTimeSpacingProduct>                               *
// *            <MCObs>...</MCObs>                                          *
// *        </ReferenceMassTimeSpacingProduct>                              *
// *        <LatticeAnisotropy>        (optional a_s/a_t) (unity if absent) *
// *            <MCObs>...</MCObs>                                          *
// *        </LatticeAnisotropy>                                            *
// *        <ParticleMass>                                                  *
// *           <Name>pion</Name> (should match names used in K-matrix)      *
// *           <MCObs>...</MCObs>                                           *
// *        </ParticleMass>                                                 *
// *           ... other particle masses                                    *
// *      </MCEnsembleParameters>                                           *
// *                                                                        *
// *    When specifying an MCObsInfo, either a short form or a long form    *
// *    can be used (must be nonsimple and real):                           *
// *                                                                        *
// *     <MCObservable>                                                     *
// *       <ObsName>T1up_Energy</ObsName> (32 char or less, no blanks)      *
// *       <Index>3</Index>        (opt nonneg integer: default 0)          *
// *     </MCObservable>                                                    *
// *                                                                        *
// *     <MCObs>T1up_Energy 3</MCObs>                                       *
// *                                                                        *
// *                                                                        *
// *    For each KB quantization block, a tag of the form is needed:        *
// *                                                                        *
// *      <KBBlock>                                                         *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <BoxQuantization>                                               *
// *          <TotalMomentumRay>ar</TotalMomentumRay>                       *
// *          <TotalMomentumIntSquared>0</TotalMomentumIntSquared>          *
// *          <LGIrrep>T1u</LGIrrep>                                        *
// *          <LmaxValues>5 3</LmaxValues>  (one for each decay channel)    *
// *        </BoxQuantization>                                              *
// *        <LabFrameEnergy>                                                *
// *           <MCObs>...</MCObs>  (record key to get data)                 *
// *        </LabFrameEnergy>                                               *
// *           .... other energies                                          *
// *      </KBBlock>                                                        *
// *                                                                        *
// *                                                                        *
// *    Specification of the K-matrix or the inverse of the K-matrix is     *
// *    done using XML of the format:                                       *
// *                                                                        *
// *      <KtildeMatrixInverse>  (or <KtildeMatrix>)                        *
// *                                                                        *
// *        <DecayChannels>                                                 *
// *           <DecayChannelInfo>                                           *
// *              <Particle1Name>pion</Particle1Name>                       *
// *              <Spin1TimesTwo>0</Spin1TimesTwo>                          *
// *              <Identical/> (if identical, do not include tags below)    *
// *              <Particle2Name>eta</Particle2Name>                        *
// *              <Spin2TimesTwo>2</Spin2TimesTwo>                          *
// *             <IntrinsicParities>same</IntrinsicParities> (or "opposite")*
// *           </DecayChannelInfo>                                          *
// *                                                                        *
// *            ... other channels infos ...                                *
// *                                                                        *
// *          (Order matters: first <DecayChannelInfo> tag is channel 0,    *
// *           second <DecayChannelInfo> is channel 1, and so on.  In the   *
// *           K-matrix, channels are referred to using the index 0, 1,...) *
// *                                                                        *
// *        </DecayChannels>                                                *
// *                                                                        *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>2</JTimesTwo>                                    *
// *            <KIndex>L(1) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(1) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *             <Polynomial><Powers>1 3</Powers></Polynomial>              *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>6</JTimesTwo>                                    *
// *            <KIndex>L(3) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(3) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *             <Polynomial><Degree>0 </Degree></Polynomial>               *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>10</JTimesTwo>                                   *
// *            <KIndex>L(5) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(5) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *            <Polynomial><Degree>0</Degree></Polynomial>                 *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *                                                                        *
// *        <StartingValues>                                                *
// *                                                                        *
// *          <KFitParamInfo>                                               *
// *            <PolynomialTerm>                                            *
// *               <Power>3</Power>                                         *
// *               <KElementInfo>...</KElementInfo>                         *
// *            </PolynomialTerm>                                           *
// *            <StartingValue>0.4534</StartingValue>                       *
// *          </KFitParamInfo>                                              *
// *          <KFitParamInfo>                                               *
// *            <PoleEnergy>                                                *
// *               <Index>3</Index>                                         *
// *               <JTimesTwo>2</JTimesTwo>                                 *
// *            </PoleEnergy>                                               *
// *            <StartingValue>2.2</StartingValue>                          *
// *          </KFitParamInfo>                                              *
// *          <KFitParamInfo>                                               *
// *            <PoleCoupling>                                              *
// *               <Index>3</Index>                                         *
// *               <JTimesTwo>2</JTimesTwo>                                 *
// *               <KIndex>...</KIndex>                                     *
// *            </PoleCoupling>                                             *
// *            <StartingValue>3.3</StartingValue>                          *
// *          </KFitParamInfo>                                              *
// *                                                                        *
// *               ....                                                     *
// *        </StartingValues>                                               *
// *                                                                        *
// *      </KtildeMatrixInverse>                                            *
// *                                                                        *
// *                                                                        *
// *                                                                        *
// *    Final note: The reference mass, the particle masses, and the        *
// *    anisotropy for each ensemble can be set to a fixed value by         *
// *    replacing the <MCObs>/<MCObservable> tag by a                       *
// *      <FixedValue>1.1123</FixedValue> tag.                              *
// *    Some groups, such as JLab, erroneously use such fixed values.       *
// *    The above feature is useful for determining how such an erroneous   *
// *    procedure effects the final error estimates on the K-matrix fix     *
// *    parameters. (NOTE: the fixed values for particle masses always      *
// *    refer to energy ratios over the reference mass.)                    *
// *                                                                        *
// *      AUXILIARY OUTPUT:                                                 *
// *                                                                        *
// *    In the constructor, if the "outfile_stub" is not empty, then the    *
// *    Ecm/mref, qcmsq/mrefsq for each channel, and the box matrix         *
// *    element samplings are written to files.  There will be one file     *
// *    for each ensemble, with suffices ".ens0", ".ens1", and so on.       *
// *                                                                        *
// *    For a given Elab energy, a particular key is used to find the       *
// *    input.  This is specified in the input XML.  Let "K" represent      *
// *    the MCObsInfo for a given Elab.  Let "Kobs" be the observable       *
// *    name of "K", and "Kindex" be the index of "K".  Then in the output  *
// *    sampling files, the keys will be                                    *
// *                                                                        *
// *      MCObsInfo("Ecm_mref["+Kobs+"]", Kindex)                           *
// *      MCObsInfo("qcmsq_mrefsq["+Kobs+"]",Kindex)                        *
// *      MCObsInfo(BoxMatString+"["+Kobs+"]",Kindex,RealPart)              *
// *      MCObsInfo(BoxMatString+"["+Kobs+"]",Kindex,ImaginaryPart)         *
// *                                                                        *
// *    where "BoxMatString" has the form                                   *
// *                                                                        *
// *      "B[momray,Psqint,IrrepB,2S,chan][2J',L',nocc'][2J,L,nocc]"        *
// *                                                                        *
// *                                                                        *
// **************************************************************************

class DeterminantResidualFit : public ChiSquare {

  KBObsHandler* KBOH;
  std::vector<BoxQuantization*> BQ;
  std::vector<RVector> Ecm_over_mref;
  std::vector<uint> ensemble_id;
  std::vector<std::vector<ComplexHermitianMatrix>> Bmat;
  KtildeMatrixCalculator* Kmat;
  KtildeInverseCalculator* Kinv;
  double omega_mu;
  std::vector<uint> nres_per_block;
  BoxQuantization::QuantCondType qctype_enum;

public:
  DeterminantResidualFit(XMLHandler& xmlin, KBObsHandler* kboh,
                         XMLHandler& xmlout,
                         const std::string& outfile_stub = "");

  virtual ~DeterminantResidualFit();

  void clear();

  virtual void guessInitialFitParamValues(std::vector<double>& fitparams) const;

  virtual void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const;

  virtual void do_output(XMLHandler& xmlout) const;

  const std::vector<KFitParamInfo>& getFitParamInfos() const;

private:
  virtual void
  evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams);

  static void read_obs(XMLHandler& xmlin, const std::string& tag, bool get_name,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset,
                       std::string& name, const MCEnsembleInfo& mcens,
                       std::map<KBObsInfo, double>& fixed_values);

  static void read_obs(XMLHandler& xmlin, const std::string& tag,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset);

  friend class TaskHandler;
};

// ************************************************************************
#endif
