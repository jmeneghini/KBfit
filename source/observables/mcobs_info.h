#ifndef MCOBS_INFO_H
#define MCOBS_INFO_H

#include "args_handler.h"
#include "xml_handler.h"

enum ComplexArg { RealPart, ImaginaryPart };

// ********************************************************************
// *                                                                  *
// *   Objects of class "MCObsInfo" store identifying information     *
// *   about one particular Monte Carlo observable.  Each observable  *
// *   must be associated with a REAL-VALUED quantity that can be     *
// *   estimated by our Monte Carlo path integrals.  These can be     *
// *   simple quantities, such as the real or imaginary part of       *
// *   a temporal correlator for one time separation, that can be     *
// *   defined on a single gauge configuration, or much more          *
// *   complicated quantities, such as a fit parameter yielding a     *
// *   stationary-state energy, determined by fitting a decaying      *
// *   exponential to a temporal correlation function.                *
// *                                                                  *
// *   An observable is termed "simple" if it can be associated with  *
// *   the integrand of a single path integral.  Simple observables   *
// *   include                                                        *
// *    (1) the real or imaginary part of a temporal correlator       *
// *                for one time separation (no vev subtraction)      *
// *    (2) the real or imaginary part of a vacuum expectation value  *
// *                of a single operator                              *
// *   Other observables are referred to as "nonsimple".              *
// *                                                                  *
// *   The class "MCObsInfo" is meant to encompass all observables    *
// *   of interest.  Observables can be classified as "primary" or    *
// *   "secondary":  "primary" refers to operator VEVs and            *
// *   correlators of field operators, which can be of type           *
// *   "BasicLapH" or "GenIrrep"; "secondary" refers to fit           *
// *   parameters, and other user-defined observables.                *
// *                                                                  *
// *   This class is meant to match that used by "SigMonD". However,  *
// *   here, we only need secondary quantities, such as energies      *
// *   and masses from fits, anisotropy, etc.                         *
// *                                                                  *
// *   "secondary" observables:                                       *
// *                                                                  *
// *     "secondary" observables are specified using a <ObsName>      *
// *     tag, unsigned integer <Index>, as well as an optional <Arg>  *
// *     tag (assumed RealPart if absent) and an optional <Simple/>   *
// *     tag (if absent, the observable is assumed to be nonsimple).  *
// *                                                                  *
// *     For these observables, there is a constructor which takes an *
// *     XMLHandler as its single argument, and another version which *
// *     takes the <ObsName> string, <Index> integer, and so on,      *
// *     Construction by XML content requires XML in the              *
// *     following format:                                            *
// *                                                                  *
// *     <MCObservable>                                               *
// *       <ObsName>T1up_Energy</ObsName> (64 char or less, no blanks)*
// *       <Index>3</Index>        (opt nonneg integer: default 0)    *
// *       <Simple/>      (optional: if simple observable)            *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *                                                                  *
// *     <MCObs>T1up_Energy 3 s re</MCObs>  (default "re" or "im")    *
// *                                        (default "n" or "s")      *
// *                                                                  *
// ********************************************************************

// ********************************************************************
// *                                                                  *
// *   Implementation notes:                                          *
// *                                                                  *
// *   All observables are encoded in a std::vector<unsigned int>.    *
// *   The first unsigned integer icode[0] contains the following     *
// *   information:                                                   *
// *                                                                  *
// *   Content of icode[0]:                                           *
// *        rightmost bit:  0 --> real part,  1 --> imaginary part    *
// *   next rightmost bit:  0 --> simple,     1 --> nonsimple         *
// *   next rightmost bit:  0 --> primary,    1 --> secondary         *
// *    remaining 29 bits:                                            *
// *        if primary:                                               *
// *            0 = Vacuum, 1 = VEV,   2 = CorrelatorAtTimeInfo       *
// *        if secondary:                                             *
// *            the unsigned integer index                            *
// *                                                                  *
// *   Here, we do not allow "primary" observables.                   *
// *   For "secondary" observables, the remaining elements of the     *
// *   icode[] vector contain the maximum 64-character "ObsName"      *
// *   string converted byte-by-byte to integers by ASCII code.       *
// *                                                                  *
// *                                                                  *
// ********************************************************************

class MCObsInfo {

  std::vector<unsigned int> icode;

public:
  MCObsInfo();

  MCObsInfo(XMLHandler& xml_in);

  MCObsInfo(const std::string& obsname, uint index = 0, bool simple = false,
            ComplexArg arg = RealPart);

  MCObsInfo(const MCObsInfo& B) : icode(B.icode) {}

  MCObsInfo& operator=(const MCObsInfo& B) {
    icode = B.icode;
    return *this;
  }

  ~MCObsInfo() {}

  void setToRealPart();

  void setToImaginaryPart();

  void resetObsIndex(uint ind); // for secondary observables

  // output functions

  bool isVacuum() const;

  bool isRealPart() const;

  bool isImaginaryPart() const;

  bool isSimple() const;

  bool isNonSimple() const;

  bool isPrimary() const;

  bool isSecondary() const;

  // routines below throw exception if inappropriate

  std::string getObsName() const;

  uint getObsIndex() const;

  std::string output(bool longform = false, int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout, bool longform = false) const; // XML output

  bool operator==(const MCObsInfo& rhs) const;

  bool operator!=(const MCObsInfo& rhs) const;

  bool operator<(const MCObsInfo& rhs) const;

  // For HDF5 support - serialize method
  std::string serialize() const;

  //  Routines below are used when MCObsInfo is a record key in
  //  an IOMap.  The IOMap class requires that every record key
  //  must occupy the same number of bytes.

  static int numints() { return max_ints; }

  void copyTo(unsigned int* buf) const;

  explicit MCObsInfo(const unsigned int* buf);

  size_t numbytes() const { return max_ints * sizeof(unsigned int); }

private:
  static const unsigned int max_ints = 24;

  void encode(const std::string& name, uint index, bool simple, ComplexArg arg);

  void set_real_part();

  void set_imag_part();

  void set_index(uint ind);

  void set_arg(ComplexArg arg);

  void assign(ArgsHandler& xt);

  void assign_from_string(const std::string& opstring);

  std::vector<std::string> split(const std::string& astr, char delimiter) const;

  std::string get_obs_name() const;

  uint get_obs_index() const;

  explicit MCObsInfo(const std::vector<unsigned int> code) : icode(code) {}

  friend class KBObsInfo;
};

// ***************************************************************

#endif
