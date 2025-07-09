#include "chisq_detres.h"
#include "chisq_spectrum.h"
#include "chisq_fit.h"
#include "task_handler.h"
#include <mpi.h>

using namespace std;

// ****************************************************************************************************
// *                                                                                                  *
// *    XML FORMAT FOR CHI-SQUARE FITTING (DoFit)                                                     *
// *                                                                                                  *
// *    This document details the complete XML structure for the <DoFit> task. This task performs a   *
// *    chi-square minimization to fit scattering parameters (via a K-matrix) and, optionally,        *
// *    lattice QCD parameters.                                                                       *
// *                                                                                                  *
// ****************************************************************************************************
// *                                                                                                  *
// *    ============================================================================================  *
// *    I. General Task Structure                                                                     *
// *    ============================================================================================  *
// *    The task requires a specific structure, starting with the action and fit type, followed by    *
// *    optional configurations and a mandatory, method-specific configuration block.                 *
// *                                                                                                  *
// *    <Task>                                                                                        *
// *      <Action>DoFit</Action>                                                                      *
// *      <Type>DeterminantResidualFit</Type> (or SpectrumFit)                                        *
// *                                                                                                  *
// *      <!-- Optional top-level configurations -->                                                  *
// *      <MinimizerInfo> ... </MinimizerInfo>                                                         *
// *      <OutSamplingsFile> ... </OutSamplingsFile>                                                   *
// *      <EcmQcmBoxSamplingsStub> ... </EcmQcmBoxSamplingsStub>                                       *
// *                                                                                                  *
// *      <!-- Mandatory method-specific configuration block -->                                      *
// *      <DeterminantResidualFit> ... </DeterminantResidualFit>                                       *
// *          or                                                                                      *
// *      <SpectrumFit> ... </SpectrumFit>                                                             *
// *    </Task>                                                                                       *
// *                                                                                                  *
// *    --------------------------------------------------------------------------------------------  *
// *    Top-Level Tag Details                                                                         *
// *    --------------------------------------------------------------------------------------------  *
// *    <Type>: (Mandatory) Specifies the fitting method.                                             *
// *      - "DeterminantResidualFit": Uses the Determinant Residual method where residuals are        *
// *        formed from det(Kinv - B) or similar quantities.                                          *
// *      - "SpectrumFit": Fits K-matrix and lattice parameters by directly comparing predicted       *
// *        and measured energy levels, including lattice parameter priors in the chi-squared.        *
// *                                                                                                  *
// *    <MinimizerInfo>: (Optional) Configures the Minuit2 minimizer.                                 *
// *      <Method>: Name of the minimization algorithm. Default: "Minuit2Migrad".                     *
// *      <ParameterRelTol>: Relative tolerance for fit parameters. Default: 1e-6.                    *
// *      <ChiSquareRelTol>: Relative tolerance for chi-squared value. Default: 1e-4.                 *
// *      <MaximumIterations>: Max number of iterations. Default: 1024.                               *
// *      <Verbosity>: Logging level. Can be "Low", "Medium", "High". Default: "Low".                 *
// *                                                                                                  *
// *    <OutSamplingsFile>: (Mandatory) Path to the output HDF5 file where the samplings of the        *
// *       best-fit parameters will be stored.                                                        *
// *                                                                                                  *
// *    <EcmQcmBoxSamplingsStub>: (Optional, for DetRes fit only) A file stub for outputting          *
// *       intermediate diagnostic quantities like Ecm/mref, qcm^2/mref^2, and Box matrix elements.    *
// *                                                                                                  *
// ****************************************************************************************************
// *                                                                                                  *
// *    ============================================================================================  *
// *    II. Fit Configuration Block (<DeterminantResidualFit> or <SpectrumFit>)                       *
// *    ============================================================================================  *
// *    This block contains the core configuration. It has method-specific tags and common tags.      *
// *                                                                                                  *
// *    --------------------------------------------------------------------------------------------  *
// *    A. Method-Specific Tags                                                                       *
// *    --------------------------------------------------------------------------------------------  *
// *    <RootFinder>: (Mandatory for SpectrumFit) Configures the numerical root-finding algorithm.    *
// *      <MaxIterations>: (Optional) Max iterations for the root finder. Default: 100.               *
// *      <InitialBracketFactor>: (Optional) Factor to expand search bracket. Default: 1.2.           *
// *      <Tolerance>: (Optional) The y-value tolerance for the root. Default: 1e-9.                  *
// *      <XTolerance>: (Optional) The x-value (energy) tolerance for the root. Default: 1e-9.        *
// *                                                                                                  *
// *    --------------------------------------------------------------------------------------------  *
// *    B. Common Tags (used inside the fit configuration block)                                      *
// *    --------------------------------------------------------------------------------------------  *
// *    <QuantizationCondition>: (Mandatory) Specifies the QC to use, e.g., "KtildeinvB", "StildeCB". *
// *                                                                                                  *
// *    <OmegaMu>: (Optional) Specifies the mu parameter for Omega-based QCs (like StildeCB).         *
// *                                                                                                  *
// *    <Verbose/>: (Optional) A flag to enable verbose logging during the fit setup.                 *
// *                                                                                                  *
// *    <DefaultEnergyFormat>: (Optional) Specifies units for energy/mass inputs. Can be              *
// *      "reference_ratio" (default) or "time_spacing_product".                                      *
// *                                                                                                  *
// *    <KtildeMatrixInverse> or <KtildeMatrix>: (Mandatory) Defines the K-matrix.                    *
// *      <DecayChannels>: Defines the scattering channels.                                           *
// *        <DecayChannelInfo>: One tag for each channel, in order (0, 1, ...).                       *
// *          <Particle1Name>: Name of the first particle.                                            *
// *          <Spin1TimesTwo>: Two times the spin of the first particle.                              *
// *          <Particle2Name>: Name of the second particle.                                           *
// *          <Spin2TimesTwo>: Two times the spin of the second particle.                             *
// *          <Identical/>: Flag if particles are identical. If present, Particle2 tags are omitted.  *
// *          <IntrinsicParities>: "same" or "opposite".                                              *
// *      <Element>: Defines a single K-matrix element.                                               *
// *        <KElementInfo>: Specifies which matrix element this is.                                   *
// *          <JTimesTwo>: Two times the total angular momentum J.                                    *
// *          <KIndex>: Two needed. Specifies row/col via "L(l) 2S(s) chan(c)" format.                *
// *        <FitForm>: The analytical form of this matrix element.                                    *
// *          <Polynomial>: A simple polynomial in Ecm.                                               *
// *            <Degree>: The highest power, e.g., "2" for c0 + c1*x + c2*x^2.                        *
// *            <Powers>: Specific powers, e.g., "0 2 4" for c0 + c1*x^2 + c2*x^4.                    *
// *      <StartingValues>: Provides initial guesses for all fit parameters in the K-matrix.          *
// *        <KFitParamInfo>: One for each parameter.                                                  *
// *          <PolynomialTerm> or <PoleEnergy> or <PoleCoupling> etc: Identifies the parameter.       *
// *          <StartingValue>: The initial floating-point value for the parameter.                    *
// *                                                                                                  *
// *    <MCEnsembleParameters>: (Mandatory) One block for each Monte Carlo ensemble.                  *
// *      <MCEnsembleInfo>: Identifies the ensemble, matching the input data file.                    *
// *      <ReferenceMassTimeSpacingProduct>: The lattice scale (a_t*m_ref). Takes:                    *
// *        <MCObs>: An observable key, e.g., "mref_times_at 0".                                      *
// *      <LatticeAnisotropy>: (Optional) The anisotropy (a_s/a_t). Takes <MCObs> or <FixedValue>.    *
// *                                      defaults to 1.0 (isotropic lattice).                        *
// *      <ParticleMass>: One for each unique decay particle. Takes <MCObs> or <FixedValue>.          *
// *        <Name>: Name of the particle (e.g., "pion").                                              *
// *                                                                                                  *
// *    <KBBlock>: (Mandatory) One for each kinematic block (defined by total momentum and irrep).    *
// *      <MCEnsembleInfo>: Associates this block with an ensemble defined above.                     *
// *      <BoxQuantization>: Defines the kinematics.                                                  *
// *        <TotalMomentumRay>: Direction of total momentum (e.g., "ar", "pd").                      *
// *        <TotalMomentumIntSquared>: Squared integer momentum, e.g., "0", "1".                      *
// *        <LGIrrep>: The little group irrep (e.g., "A1", "T1u").                                    *
// *        <LmaxValues>: Max orbital angular momentum values, one for each channel.                  *
// *      <!-- Energy/Residual Input (Choose one based on fit type) -->                               *
// *      <LabFrameEnergy>: (For DetRes Fit) The measured lab-frame energy level. Takes <MCObs>.      *
// *      <LabFrameEnergyShift>: (For Spectrum Fit) The measured energy shift. Takes <MCObs>.         *
// *      <LabFrameEnergyMin> and <LabFrameEnergyMax>: (For Spectrum Fit) Mandatory search window.    *
// *                                                                                                  *
// *    <KBObservables>: (Mandatory) Provides information about the input data files.                 *
// *      <MCSamplingInfo>: Specifies the resampling method.                                          *
// *        <Jackknife/> or <Bootstrapper><NumberResamplings>...</Bootstrapper>                       *
// *      <SamplingData>: Contains the list of input data files.                                      *
// *        <FileName>path/to/data.hdf5</FileName>                                                    *
// *                                                                                                  *
// ****************************************************************************************************

void TaskHandler::doFit(XMLHandler& xmltask, XMLHandler& xmlout,
                        int taskcount) {
  ChiSquareMinimizerInfo mz_info; // default minimizer info
  if (xmltask.count_among_children("MinimizerInfo") > 0) {
    ChiSquareMinimizerInfo mz_user(xmltask);
    mz_info = mz_user;
  }
  string outsampfile;
  xmlread(xmltask, "OutSamplingsFile", outsampfile, "doFit");
  string EcmQcmBoxSampStub;
  xmlreadif(xmltask, "EcmQcmBoxSamplingsStub", EcmQcmBoxSampStub, "doFit");
  string fittype;
  xmlreadchild(xmltask, "Type", fittype, "DoFit");
  xmlout.set_root("DoFit");
  XMLHandler xmlmz;
  mz_info.output(xmlmz);
  xmlout.put_child(xmlmz);
  xmlout.put_child("Type", fittype);
  xmlout.put_child("OutSamplingsFile", outsampfile);
  if (!EcmQcmBoxSampStub.empty())
    xmlout.put_child("EcmQcmBoxSamplingsStub", EcmQcmBoxSampStub);
  double chisq_dof, qual;
  vector<MCEstimate> bestfit_params;
  RealSymmetricMatrix param_covariance;

  if (fittype == "DeterminantResidualFit") {
    try {
      XMLHandler xmlf(xmltask, "DeterminantResidualFit");

      // get the quantization condition for the output directory
      string qctype;
      xmlreadif(xmltask, "QuantizationCondition", qctype, "DeterminantResidualFit");
      filesystem::path output_path = createKBOutputDirectory(m_output_directory, qctype);

      EcmQcmBoxSampStub = (output_path / EcmQcmBoxSampStub).string();
      outsampfile = (output_path / outsampfile).string();

      XMLHandler xmlcon;
      DeterminantResidualFit DRF(xmlf, m_obs, xmlcon, EcmQcmBoxSampStub);
      xmlout.put_child(xmlcon);
      XMLHandler xmlof;
      DRF.do_output(xmlof);
      xmlout.put_child(xmlof);
      doChiSquareFitting(DRF, mz_info, chisq_dof, qual, bestfit_params,
                         param_covariance, outsampfile, xmlout, m_obs, MPI_COMM_WORLD);
    } catch (const std::exception& xp) {
      string msg("DetRes fit failed: ");
      msg += xp.what();
      xmlout.put_child("Error", msg);
      throw(std::invalid_argument(msg));
    }
  }

  else if (fittype == "SpectrumFit") {
    try {
      XMLHandler xmlf(xmltask, "SpectrumFit");

      // get the quantization condition for the output directory
      string qctype;
      xmlreadif(xmltask, "QuantizationCondition", qctype, "SpectrumFit");
      filesystem::path output_path = createKBOutputDirectory(m_output_directory, qctype);

      EcmQcmBoxSampStub = (output_path / EcmQcmBoxSampStub).string();
      outsampfile = (output_path / outsampfile).string();

      XMLHandler xmlcon;
      SpectrumFit SF(xmlf, m_obs, xmlcon, EcmQcmBoxSampStub);
      xmlout.put_child(xmlcon);
      XMLHandler xmlof;
      SF.do_output(xmlof);
      xmlout.put_child(xmlof);
      doChiSquareFitting(SF, mz_info, chisq_dof, qual, bestfit_params,
                         param_covariance, outsampfile, xmlout, m_obs, MPI_COMM_WORLD);
    } catch (const std::exception& xp) {
      string msg("Spectrum fit failed: ");
      msg += xp.what();
      xmlout.put_child("Error", msg);
      throw(std::invalid_argument(msg));
    }
  }
}

// ***************************************************************************************
