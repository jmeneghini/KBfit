#include "box_quant.h"
#include "utils.h"

using namespace std;

// ****************************************************************

BoxQuantBasisState::BoxQuantBasisState(BoxMatrix* boxmat, uint channel_index,
                                       uint Stimestwo, uint Jtimestwo, uint L,
                                       uint occurrence)
    : m_boxmat(boxmat) {
  check_encode(channel_index, Stimestwo, Jtimestwo, L, occurrence);
}

//   m_store encodes the channel index (6 bits),
//   total S times two (6 bits), total J times two (7 bits)
//   L (7 bits), occurrence index (6 bits)

void BoxQuantBasisState::check_encode(uint channel_index, uint Stimestwo,
                                      uint Jtimestwo, uint L, uint occurrence) {
  if ((channel_index >= 64) || (Stimestwo >= 64) || (Jtimestwo >= 128) ||
      (L >= 128) || (occurrence >= 64))
    throw(std::invalid_argument("Unsupported BoxQuantBasisState value"));
  m_store = channel_index;
  m_store <<= 6;
  m_store |= Stimestwo;
  m_store <<= 7;
  m_store |= Jtimestwo;
  m_store <<= 7;
  m_store |= L;
  m_store <<= 6;
  m_store |= occurrence;
}

string BoxQuantBasisState::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string BoxQuantBasisState::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void BoxQuantBasisState::output(XMLHandler& xmlout) const {
  xmlout.set_root("BoxQuantBasisState");
  string ststr(" chan(");
  ststr += make_string(getChannelIndex()) + ")";
  ststr += " 2S(" + make_string(getStimestwo()) + ")";
  ststr += " 2J(" + make_string(getJtimestwo()) + ")";
  ststr += " L(" + make_string(getL()) + ")";
  ststr += " n(" + make_string(getOccurrence()) + ") ";
  xmlout.put_child("Ket", ststr);
}

// ****************************************************************

// *     <BoxQuantization>
// *       <TotalMomentumRay>oa</TotalMomentumRay>
// *           ["ar"=at rest (0,0,0), "oa"=on axis (0,0,1),
// *            "pd"=planar diagonal (0,1,1),
// *            "cd"=cubic diagonal (1,1,1)]
// *       <TotalMomentumIntSquared>2</TotalMomentumIntSquared>
// *           [if total momentum is (2*Pi/L)(1,1,1), then IntSquared = 3]
// *       <LGIrrep>T1u</LGIrrep>
// *         <LmaxValues>4 4</LmaxValues>
// *     </BoxQuantization>
// *

BoxQuantization::BoxQuantization(XMLHandler& xmlin,
                                 KtildeMatrixCalculator* Kmatptr)
    : m_Kmat(Kmatptr), m_Kinv(0) {
  xmlinitialize(xmlin);
}

BoxQuantization::BoxQuantization(XMLHandler& xmlin,
                                 KtildeInverseCalculator* Kinvptr)
    : m_Kmat(0), m_Kinv(Kinvptr) {
  xmlinitialize(xmlin);
}

BoxQuantization::BoxQuantization(XMLHandler& xmlin,
                                 KtildeMatrixCalculator* Kmatptr,
                                 KtildeInverseCalculator* Kinvptr)
    : m_Kmat(Kmatptr), m_Kinv(Kinvptr) {
  if ((m_Kmat != 0) && (m_Kinv != 0))
    throw(std::invalid_argument(
        string("One of the pointers must be null in BoxQuantization "
               "constructor with ") +
        string("KtildeMatrixCalculator and KtildeInverseCalculator pointers")));
  xmlinitialize(xmlin);
}

void BoxQuantization::xmlinitialize(XMLHandler& xmlin) {
  try {
    XMLHandler xmlinfo(xmlin, "BoxQuantization");
    xmlreadchild(xmlinfo, "TotalMomentumRay", m_momray);
    uint mom_int_sq = 0;
    xmlreadchild(xmlinfo, "TotalMomentumIntSquared", mom_int_sq);
    set_dvector(mom_int_sq);
    xmlreadchild(xmlinfo, "LGIrrep", m_lgirrep);
    uint nchan = getNumberOfDecayChannels();
    xmlreadchild(xmlinfo, "LmaxValues", m_Lmaxes);
    if (m_Lmaxes.size() != nchan) {
      throw(std::invalid_argument(
          "Size mismatch between decay infos and L-maxes"));
    }
    initialize();
  } catch (const std::exception& xp) {
    clear();
    throw(std::invalid_argument(string("Invalid input to BoxQuantization: ") +
                                xp.what()));
  }
}

BoxQuantization::BoxQuantization(
    const std::string& mom_ray, uint mom_int_sq, const std::string& lgirrep,
    const std::vector<DecayChannelInfo>& chan_infos,
    const std::vector<uint> lmaxes, KtildeMatrixCalculator* Kmatptr)
    : m_lgirrep(lgirrep), m_momray(mom_ray), m_Lmaxes(lmaxes), m_Kmat(Kmatptr),
      m_Kinv(0) {
  set_dvector(mom_int_sq);
  if (Kmatptr->getNumberOfDecayChannels() != m_Lmaxes.size())
    throw(
        std::invalid_argument("Size mismatch between decay infos and L-maxes"));
  initialize();
}

BoxQuantization::BoxQuantization(
    const std::string& mom_ray, uint mom_int_sq, const std::string& lgirrep,
    const std::vector<DecayChannelInfo>& chan_infos,
    const std::vector<uint> lmaxes, KtildeInverseCalculator* Kinvptr)
    : m_lgirrep(lgirrep), m_momray(mom_ray), m_Lmaxes(lmaxes), m_Kmat(0),
      m_Kinv(Kinvptr) {
  set_dvector(mom_int_sq);
  if (Kinvptr->getNumberOfDecayChannels() != m_Lmaxes.size())
    throw(
        std::invalid_argument("Size mismatch between decay infos and L-maxes"));
  initialize();
}

void BoxQuantization::set_dvector(uint mom_int_sq) {
  if (m_momray == "ar") { // at rest
    if (mom_int_sq != 0)
      throw(
          std::invalid_argument("At rest ray with nonzero momentum squared!"));
    m_dx = m_dy = m_dz = 0;
  } else if (m_momray == "oa") { // on-axis
    m_dx = m_dy = 0;
    if (mom_int_sq == 1)
      m_dz = 1;
    else if (mom_int_sq == 4)
      m_dz = 2;
    else if (mom_int_sq == 9)
      m_dz = 3;
    else if (mom_int_sq == 16)
      m_dz = 4;
    else if (mom_int_sq == 25)
      m_dz = 5;
    else if (mom_int_sq == 36)
      m_dz = 6;
    else
      throw(std::invalid_argument("Supported momentum squared"));
  } else if (m_momray == "pd") { // planar-diagonal
    if (mom_int_sq == 2)
      m_dz = 1;
    else if (mom_int_sq == 8)
      m_dz = 2;
    else if (mom_int_sq == 18)
      m_dz = 3;
    else if (mom_int_sq == 32)
      m_dz = 4;
    else if (mom_int_sq == 50)
      m_dz = 5;
    else if (mom_int_sq == 72)
      m_dz = 6;
    else
      throw(std::invalid_argument("Supported momentum squared"));
    m_dx = 0;
    m_dy = m_dz;
  } else if (m_momray == "cd") { // cubic-diagonal
    if (mom_int_sq == 3)
      m_dz = 1;
    else if (mom_int_sq == 12)
      m_dz = 2;
    else if (mom_int_sq == 27)
      m_dz = 3;
    else if (mom_int_sq == 48)
      m_dz = 4;
    else if (mom_int_sq == 75)
      m_dz = 5;
    else if (mom_int_sq == 108)
      m_dz = 6;
    else
      throw(std::invalid_argument("Supported momentum squared"));
    m_dx = m_dy = m_dz;
  } else {
    throw(std::invalid_argument("Unsupported momentum ray"));
  }
}

void BoxQuantization::clear() {
  m_Lmaxes.clear();
  m_masses1.clear();
  m_masses2.clear();
  for (list<pair<BoxMatrix*, uint>>::iterator it = m_boxes.begin();
       it != m_boxes.end(); it++)
    delete it->first;
  m_boxes.clear();
  for (list<WZetaRGLCalculator*>::iterator it = m_wzetas.begin();
       it != m_wzetas.end(); it++)
    delete (*it);
  m_wzetas.clear();
  m_basis.clear();
}

BoxQuantization::~BoxQuantization() { clear(); }

BoxQuantization::BoxQuantization() 
    : m_mref_L(6.0), m_Kmat(nullptr), m_Kinv(nullptr) {
  // Initialize to default/empty state - will be populated by clone method
}

std::unique_ptr<BoxQuantization> BoxQuantization::clone(KtildeMatrixCalculator* new_Kmat,
                                                       KtildeInverseCalculator* new_Kinv) const {
  // Create new instance using default constructor
  auto cloned = std::unique_ptr<BoxQuantization>(new BoxQuantization());
  
  // Copy basic value-type members
  cloned->m_lgirrep = this->m_lgirrep;
  cloned->m_lgirrepB = this->m_lgirrepB;
  cloned->m_momray = this->m_momray;
  cloned->m_dx = this->m_dx;
  cloned->m_dy = this->m_dy;
  cloned->m_dz = this->m_dz;
  cloned->m_Lmaxes = this->m_Lmaxes;
  cloned->m_mref_L = this->m_mref_L;
  cloned->m_masses1 = this->m_masses1;
  cloned->m_masses2 = this->m_masses2;
  cloned->m_decay_channel_S_num = this->m_decay_channel_S_num;
  
  // Copy basis (value types in set)
  cloned->m_basis = this->m_basis;
  
  // Set K-matrix calculator pointers (external objects - shallow copy unless new ones provided)
  if (new_Kmat != nullptr) {
    cloned->m_Kmat = new_Kmat;
    cloned->m_Kinv = nullptr;
  } else if (new_Kinv != nullptr) {
    cloned->m_Kinv = new_Kinv;
    cloned->m_Kmat = nullptr;
  } else {
    // Use existing pointers (shallow copy)
    cloned->m_Kmat = this->m_Kmat;
    cloned->m_Kinv = this->m_Kinv;
  }
  
  // Deep copy WZetaRGLCalculator objects
  for (const auto* wzeta : this->m_wzetas) {
    if (wzeta != nullptr) {
      // WZetaRGLCalculator has private copy constructor, so create new one with same parameters
      WZetaRGLCalculator* cloned_wzeta = new WZetaRGLCalculator();
      cloned_wzeta->reset(wzeta->getSVector(), wzeta->getGamma(), wzeta->getUsq());
      cloned->m_wzetas.push_back(cloned_wzeta);
    } else {
      cloned->m_wzetas.push_back(nullptr);
    }
  }
  
  // Deep copy BoxMatrix objects
  // Note: We need to update the BoxQuantBasisState objects to point to the new BoxMatrix objects
  std::map<BoxMatrix*, BoxMatrix*> box_ptr_map; // old -> new mapping
  
  for (const auto& pair : this->m_boxes) {
    BoxMatrix* original_box = pair.first;
    uint channel_index = pair.second;
    
    if (original_box != nullptr) {
      // BoxMatrix has private copy constructor, so create new one with same parameters
      // Need to find the corresponding WZetaRGLCalculator for this box
      // The order of boxes and wzetas should correspond based on how they were created
      auto wzeta_it = cloned->m_wzetas.begin();
      std::advance(wzeta_it, std::distance(this->m_boxes.begin(), 
                                          std::find_if(this->m_boxes.begin(), this->m_boxes.end(),
                                                      [original_box](const auto& p) { return p.first == original_box; })));
      
      BoxMatrix* cloned_box = new BoxMatrix(original_box->getEcmTransform(), 
                                          **wzeta_it,
                                          original_box->getTotalSpinTimesTwo(),
                                          original_box->getLittleGroupIrrep(),
                                          original_box->getLmax());
      cloned->m_boxes.push_back(std::make_pair(cloned_box, channel_index));
      box_ptr_map[original_box] = cloned_box;
    } else {
      cloned->m_boxes.push_back(std::make_pair(nullptr, channel_index));
      box_ptr_map[original_box] = nullptr;
    }
  }
  
  // Update basis states to point to new BoxMatrix objects
  cloned->m_basis.clear();
  for (const auto& basis_state : this->m_basis) {
    BoxMatrix* old_ptr = basis_state.getBoxMatrixPtr();
    BoxMatrix* new_ptr = box_ptr_map[old_ptr];
    
    BoxQuantBasisState new_state(new_ptr, 
                                 basis_state.getChannelIndex(),
                                 basis_state.getStimestwo(),
                                 basis_state.getJtimestwo(),
                                 basis_state.getL(),
                                 basis_state.getOccurrence());
    cloned->m_basis.insert(new_state);
  }
  
  return cloned;
}

void BoxQuantization::initialize() {
  // set up default values of masses; these can be changed later
  const vector<DecayChannelInfo>& decay_infos = getDecayChannelInfos();
  uint nchan = decay_infos.size();
  if (nchan == 0)
    throw(std::invalid_argument("Must have at least ONE decay channel"));
  bool ipprod = decay_infos[0].sameIntrinsicParities();
  for (uint a = 1; a < nchan; a++) {
    if (ipprod != decay_infos[a].sameIntrinsicParities())
      throw(std::invalid_argument(
          "All channels must have same product of intrinsic parities"));
  }
  if (ipprod)
    m_lgirrepB = m_lgirrep;
  else
    m_lgirrepB = getLGIrrepParityFlip();
  m_mref_L = 6.0; // arbitrary default value; should be reset
  m_masses1.resize(nchan);
  m_masses2.resize(nchan);
  for (uint k = 0; k < nchan; k++)
    m_masses1[k] = m_masses2[k] = 1.0; // arbitrary default values
  // set up the basis states
  setup_basis();
}

std::optional<BoxQuantization::QuantCondType>
BoxQuantization::getQuantCondTypeFromString(const std::string& qctype) const {
  // convert string to lower case
  std::string qctype_lower = qctype;
  std::transform(qctype_lower.begin(), qctype_lower.end(), qctype_lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  auto it = m_qctype_map.find(qctype_lower);
  return (it == m_qctype_map.end()) ? std::nullopt
                                    : std::make_optional(it->second);
}

void BoxQuantization::setRefMassL(double mref_L) {
  if (mref_L <= 0.0)
    throw(
        std::invalid_argument("Unsupported reference mass: zero or negative"));
  m_mref_L = mref_L;
  const vector<DecayChannelInfo>& decay_infos = getDecayChannelInfos();
  for (list<pair<BoxMatrix*, uint>>::iterator it = m_boxes.begin();
       it != m_boxes.end(); it++) {
    uint chan = it->second;
    if (decay_infos[chan].areIdentical())
      it->first->resetMasses(m_mref_L, m_masses1[chan]);
    else
      it->first->resetMasses(m_mref_L, m_masses1[chan], m_masses2[chan]);
  }
}

void BoxQuantization::setMassesOverRef(uint channel_index,
                                       double mass1_over_ref,
                                       double mass2_over_ref) {
  if ((mass1_over_ref <= 0.0) || (mass2_over_ref <= 0.0))
    throw(
        std::invalid_argument("Unsupported particle masses: zero or negative"));
  m_masses1.at(channel_index) = mass1_over_ref;
  m_masses2.at(channel_index) = mass2_over_ref;
  const vector<DecayChannelInfo>& decay_infos = getDecayChannelInfos();
  for (list<pair<BoxMatrix*, uint>>::iterator it = m_boxes.begin();
       it != m_boxes.end(); it++) {
    uint chan = it->second;
    if (chan == channel_index) {
      if (decay_infos[chan].areIdentical()) {
        if (mass1_over_ref != mass2_over_ref)
          throw(std::invalid_argument(
              "Identical particles cannot have different masses"));
        it->first->resetMasses(m_mref_L, m_masses1[chan]);
      } else
        it->first->resetMasses(m_mref_L, m_masses1[chan], m_masses2[chan]);
    }
  }
}

const EcmTransform&
BoxQuantization::getDecayChannelEcmTransform(uint channel_index) const {
  if (channel_index >= m_decay_channel_S_num.size())
    throw std::out_of_range("channel_index too large");

  auto box_it = m_boxes.begin();
  std::size_t skip = 0;
  for (std::size_t c = 0; c < channel_index; ++c)
    skip += m_decay_channel_S_num[c];

  std::advance(box_it, skip);
  return box_it->first->getEcmTransform();
}

string BoxQuantization::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string BoxQuantization::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void BoxQuantization::output(XMLHandler& xmlout) const {
  xmlout.set_root("BoxQuantization");
  xmlout.put_child("TotalMomentumRay", getMomRay());
  xmlout.put_child("TotalMomentumIntSquared",
                   make_string(getTotalMomentumIntegerSquared()));
  xmlout.put_child("LGIrrep", getLittleGroupIrrep());
  xmlout.put_child("LmaxValues", make_string(m_Lmaxes));
}

void BoxQuantization::outputBasis(XMLHandler& xmlout) const {
  xmlout.set_root("BoxQuantizationBasis");
  for (std::set<BoxQuantBasisState>::iterator it = m_basis.begin();
       it != m_basis.end(); it++) {
    XMLHandler xmlb;
    it->output(xmlb);
    xmlout.put_child(xmlb);
  }
}

string BoxQuantization::outputBasis(int indent) const {
  XMLHandler xmlout;
  outputBasis(xmlout);
  return xmlout.output(indent);
}

// returns key with name "B[momray,Psqint,IrrepB,2S,chan][2J',L',n'][2J,L,n]"

std::string BoxQuantization::getKeyString(int row, int col) const {
  if ((row < 0) || (col < 0) || (row >= int(m_basis.size())) ||
      (col >= int(m_basis.size())))
    throw(std::invalid_argument("Bad indices in BoxQuantization::getKey"));
  std::set<BoxQuantBasisState>::iterator itrow, itcol;
  itrow = m_basis.begin();
  itcol = m_basis.begin();
  for (int k = 0; k < row; ++k)
    ++itrow;
  for (int k = 0; k < col; ++k)
    ++itcol;
  uint chanrow = itrow->getChannelIndex();
  uint twoSrow = itrow->getStimestwo();
  uint twoJrow = itrow->getJtimestwo();
  uint Lrow = itrow->getL();
  uint noccrow = itrow->getOccurrence();
  uint chancol = itcol->getChannelIndex();
  uint twoScol = itcol->getStimestwo();
  uint twoJcol = itcol->getJtimestwo();
  uint Lcol = itcol->getL();
  uint nocccol = itcol->getOccurrence();
  if ((chanrow != chancol) || (twoSrow != twoScol))
    throw(std::invalid_argument("Bad S,a in BoxQuantization::getKey"));
  string keyname("B[");
  keyname += getMomRay() + "," + make_string(getTotalMomentumIntegerSquared()) +
             "," + getLittleGroupBoxIrrep();
  keyname += "," + make_string(twoSrow) + "," + make_string(chanrow) + "]";
  keyname += "[" + make_string(twoJrow) + "," + make_string(Lrow) + "," +
             make_string(noccrow) + "]";
  keyname += "[" + make_string(twoJcol) + "," + make_string(Lcol) + "," +
             make_string(nocccol) + "]";
  return keyname;
}

void BoxQuantization::outputKFitParams(XMLHandler& xmlout) const // XML output
{
  xmlout.set_root("BoxQuantizationKFitParams");
  const std::vector<KFitParamInfo>& fref = getFitParameterInfos();
  for (std::vector<KFitParamInfo>::const_iterator it = fref.begin();
       it != fref.end(); it++) {
    XMLHandler xmlk;
    it->output(xmlk);
    xmlout.put_child(xmlk);
  }
}

std::string BoxQuantization::outputKFitParams(int indent) const // XML output
{
  XMLHandler xmlout;
  outputKFitParams(xmlout);
  return xmlout.output(indent);
}

const std::vector<KFitParamInfo>&
BoxQuantization::getFitParameterInfos() const {
  if (m_Kmat)
    return m_Kmat->getFitParameterInfos();
  else
    return m_Kinv->getFitParameterInfos();
}

int BoxQuantization::getParameterIndex(
    const KFitParamInfo& kinfo) const // returns -1 if not found
{
  if (m_Kmat)
    return m_Kmat->getParameterIndex(kinfo);
  else
    return m_Kinv->getParameterIndex(kinfo);
}

double BoxQuantization::getParameterValue(const KFitParamInfo& kinfo) const {
  if (m_Kmat)
    return m_Kmat->getParameterValue(kinfo);
  else
    return m_Kinv->getParameterValue(kinfo);
}

set<KElementInfo> BoxQuantization::getElementInfos() const {
  if (m_Kmat)
    return m_Kmat->getElementInfos();
  else
    return m_Kinv->getElementInfos();
}

string BoxQuantization::getLGIrrepParityFlip() {
  if (m_momray == "ar") {
    string res(m_lgirrep);
    uint n = res.length();
    if (res[n - 1] == 'g')
      res[n - 1] = 'u';
    else if (res[n - 1] == 'u')
      res[n - 1] = 'g';
    else
      throw(std::invalid_argument("Unsupported irrep"));
    return res;
  } else if (m_momray == "oa") {
    if (m_lgirrep == "A1")
      return string("A2");
    else if (m_lgirrep == "A2")
      return string("A1");
    else if (m_lgirrep == "B1")
      return string("B2");
    else if (m_lgirrep == "B2")
      return string("B1");
    else if ((m_lgirrep == "E") || (m_lgirrep == "G1") || (m_lgirrep == "G2"))
      return m_lgirrep;
    else
      throw(std::invalid_argument("Unsupported irrep"));
  } else if (m_momray == "pd") {
    if (m_lgirrep == "A1")
      return string("A2");
    else if (m_lgirrep == "A2")
      return string("A1");
    else if (m_lgirrep == "B1")
      return string("B2");
    else if (m_lgirrep == "B2")
      return string("B1");
    else if (m_lgirrep == "G")
      return m_lgirrep;
    else
      throw(std::invalid_argument("Unsupported irrep"));
  } else if (m_momray == "cd") {
    if (m_lgirrep == "A1")
      return string("A2");
    else if (m_lgirrep == "A2")
      return string("A1");
    else if (m_lgirrep == "F1")
      return string("F2");
    else if (m_lgirrep == "F2")
      return string("F1");
    else if ((m_lgirrep == "E") || (m_lgirrep == "G"))
      return m_lgirrep;
    else
      throw(std::invalid_argument("Unsupported irrep"));
  } else
    throw(std::invalid_argument("Unsupported momentum ray"));
}

// *  Loop over the channels
// *  For each channel,
// *      - create a Ecmtransform object
// *      - create a WZetaRGLCalculator object
// *      - loop over allowed values of S  (|s1-s2|,...,s1+s2)
// *           - create box matrix for each
// *           - each box matrix creates the n J L basis states

// *   Basis states in the current block are labelled by
// *      irrep occurrence  J L S channel_index
// *  (irrep_row=1 always taken here)

// *   Lastly, looks for states whose K-matrix elements are
// *   all zero.  These states are excluded from the basis

void BoxQuantization::setup_basis() {
  const vector<DecayChannelInfo>& decay_infos = getDecayChannelInfos();
  uint nchan = decay_infos.size();
  for (uint chan = 0; chan < nchan; chan++) {
    EcmTransform Ecm(m_dx, m_dy, m_dz, m_mref_L, m_masses1[chan],
                     m_masses2[chan]);
    WZetaRGLCalculator* mzptr = new WZetaRGLCalculator;
    m_wzetas.push_back(mzptr);
    uint s1timestwo = decay_infos[chan].getSpin1timestwo();
    uint s2timestwo = decay_infos[chan].getSpin2timestwo();
    uint Stimestwomax =
        min(s1timestwo + s2timestwo, BoxMatrix::getTotalSpinTimesTwoMax(Ecm));
    uint Stimestwomin = std::abs(int(s1timestwo) - int(s2timestwo));
    uint num_S = (Stimestwomax - Stimestwomin)/2 + 1;
    m_decay_channel_S_num.push_back(num_S);
    for (uint Stimestwo = Stimestwomin; Stimestwo <= Stimestwomax;
         Stimestwo += 2) {
      BoxMatrix* mbptr =
          new BoxMatrix(Ecm, *mzptr, Stimestwo, m_lgirrepB, m_Lmaxes[chan]);
      uint nelem = mbptr->getNumberOfIndepElements();
      if (nelem > 0) {
        m_boxes.push_back(make_pair(mbptr, chan));
        for (uint i = 0; i < nelem; i++) {
          BoxMatrixQuantumNumbers bqm(mbptr->getQuantumNumbers(i));
          m_basis.insert(BoxQuantBasisState(mbptr, chan, Stimestwo,
                                            bqm.getRowJtimestwo(),
                                            bqm.getRowL(), bqm.getRowNocc()));
          m_basis.insert(BoxQuantBasisState(
              mbptr, chan, Stimestwo, bqm.getColumnJtimestwo(),
              bqm.getColumnL(), bqm.getColumnNocc()));
        }
      } else {
        delete mbptr;
      }
    }
  }
  if (m_basis.empty()) {
    throw(std::invalid_argument(string(
        "Null basis in BoxQuantization even before Kmatrix exclusions")));
  }
  // look for states to exclude due to K-matrix
  set<BoxQuantBasisState> exclusions;
  if (m_Kmat != 0)
    exclusions = find_excluded_states_from_ktilde(m_Kmat);
  else
    exclusions = find_excluded_states_from_ktilde(m_Kinv);
  // if ALL states excluded, this means Kmatrix is zero; end user needs to know
  // which K matrix elements should be set; throw this information
  if (exclusions.size() == m_basis.size()) {
    set<KElementInfo> kelems;
    for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
         rt != m_basis.end(); rt++)
      for (set<BoxQuantBasisState>::iterator ct = rt; ct != m_basis.end();
           ct++) {
        try {
          KElementInfo kel(rt->getJtimestwo(), rt->getL(), rt->getStimestwo(),
                           rt->getChannelIndex(), ct->getL(),
                           ct->getStimestwo(), ct->getChannelIndex());
          kelems.insert(kel);
        } catch (const std::exception& xp) {
        }
      }
    string kstr;
    for (std::set<KElementInfo>::const_iterator it = kelems.begin();
         it != kelems.end(); it++) {
      XMLHandler xmle("Element");
      XMLHandler xmli;
      it->output(xmli);
      xmle.put_child(xmli);
      xmle.put_child("FitForm", "...");
      kstr += xmle.output(1);
    }
    throw(std::invalid_argument(
        string("Null basis in BoxQuantization after Kmatrix exclusions: need K "
               "elements follow: ") +
        kstr));
  }
  for (set<BoxQuantBasisState>::iterator it = exclusions.begin();
       it != exclusions.end(); it++) {
    m_basis.erase(*it);
  }
}

double BoxQuantization::getEcmOverMrefFromElab(double Elab_over_mref) const {
  return m_boxes.front().first->getEcmOverMrefFromElab(Elab_over_mref);
}

double BoxQuantization::getElabOverMrefFromEcm(double Ecm_over_mref) const {
  return m_boxes.front().first->getElabOverMrefFromEcm(Ecm_over_mref);
}

void BoxQuantization::getQcmsqOverMrefsqFromElab(
    double Elab_over_mref, RVector& qcmsq_over_mrefsq) const {
  uint nchan = getNumberOfDecayChannels();
  qcmsq_over_mrefsq.resize(nchan);
  for (std::list<std::pair<BoxMatrix*, uint>>::const_iterator it =
           m_boxes.begin();
       it != m_boxes.end(); ++it) {
    qcmsq_over_mrefsq[it->second] =
        it->first->getQcmsqOverMrefsqFromElab(Elab_over_mref);
  }
}

void BoxQuantization::getQcmsqOverMrefsqFromEcm(
    double Ecm_over_mref, RVector& qcmsq_over_mrefsq) const {
  uint nchan = getNumberOfDecayChannels();
  qcmsq_over_mrefsq.resize(nchan);
  for (std::list<std::pair<BoxMatrix*, uint>>::const_iterator it =
           m_boxes.begin();
       it != m_boxes.end(); ++it) {
    qcmsq_over_mrefsq[it->second] =
        it->first->getQcmsqOverMrefsqFromEcm(Ecm_over_mref);
  }
}

void BoxQuantization::getBoxMatrixFromElab(double Elab_over_mref,
                                           ComplexHermitianMatrix& B) {
  CMatrix dummy;
  get_box_matrix(Elab_over_mref, B, dummy, true, true);
}

void BoxQuantization::getBoxMatrixFromEcm(double Ecm_over_mref,
                                          ComplexHermitianMatrix& B) {
  CMatrix dummy;
  get_box_matrix(Ecm_over_mref, B, dummy, false, true);
}

void BoxQuantization::getKtildeFromElab(double Elab_over_mref,
                                        RealSymmetricMatrix& Ktilde) {
  if (m_Kmat == 0)
    throw(std::runtime_error("Cannot evaluate Ktilde in Ktildeinverse mode"));
  RMatrix dummy;
  get_ktilde_matrix(Elab_over_mref, Ktilde, dummy, true, true, m_Kmat);
}

void BoxQuantization::getKtildeFromEcm(double Ecm_over_mref,
                                       RealSymmetricMatrix& Ktilde) {
  if (m_Kmat == 0)
    throw(std::runtime_error("Cannot evaluate Ktilde in Ktildeinverse mode"));
  RMatrix dummy;
  get_ktilde_matrix(Ecm_over_mref, Ktilde, dummy, false, true, m_Kmat);
}

void BoxQuantization::getKtildeinvFromElab(double Elab_over_mref,
                                           RealSymmetricMatrix& Ktildeinv) {
  if (m_Kinv == 0)
    throw(std::runtime_error("Cannot evaluate Ktildeinverse in Ktilde mode"));
  RMatrix dummy;
  get_ktilde_matrix(Elab_over_mref, Ktildeinv, dummy, true, true, m_Kinv);
}

void BoxQuantization::getKtildeinvFromEcm(double Ecm_over_mref,
                                          RealSymmetricMatrix& Ktildeinv) {
  if (m_Kinv == 0)
    throw(std::runtime_error("Cannot evaluate Ktildeinverse in Ktilde mode"));
  RMatrix dummy;
  get_ktilde_matrix(Ecm_over_mref, Ktildeinv, dummy, false, true, m_Kinv);
}

void BoxQuantization::getKtildeOrInverseFromElab(
    double Elab_over_mref, RealSymmetricMatrix& KtildeOrInverse) {
  RMatrix dummy;
  if (m_Kmat != 0)
    get_ktilde_matrix(Elab_over_mref, KtildeOrInverse, dummy, true, true,
                      m_Kmat);
  else
    get_ktilde_matrix(Elab_over_mref, KtildeOrInverse, dummy, true, true,
                      m_Kinv);
}

void BoxQuantization::getKtildeOrInverseFromEcm(
    double Ecm_over_mref, RealSymmetricMatrix& KtildeOrInverse) {
  RMatrix dummy;
  if (m_Kmat != 0)
    get_ktilde_matrix(Ecm_over_mref, KtildeOrInverse, dummy, false, true,
                      m_Kmat);
  else
    get_ktilde_matrix(Ecm_over_mref, KtildeOrInverse, dummy, false, true,
                      m_Kinv);
}

void BoxQuantization::outputKBMatricesFromElab(double Elab_over_mref,
                                               std::ostream& fout) {
  output_matrices(Elab_over_mref, true, fout);
}

void BoxQuantization::outputKBMatricesFromEcm(double Ecm_over_mref,
                                              std::ostream& fout) {
  output_matrices(Ecm_over_mref, false, fout);
}

// gets the Cayley transformed (CT) matrices

void BoxQuantization::getCTBoxMatrixFromElab(double Elab_over_mref,
                                             CMatrix& CB) {
  ComplexHermitianMatrix B;
  getBoxMatrixFromElab(Elab_over_mref, B);
  calcCayleyTransformMatrix(B, CB, true);
}

void BoxQuantization::getCTBoxMatrixFromEcm(double Ecm_over_mref, CMatrix& CB) {
  ComplexHermitianMatrix B;
  getBoxMatrixFromEcm(Ecm_over_mref, B);
  calcCayleyTransformMatrix(B, CB, true);
}

void BoxQuantization::getStildeFromElab(double Elab_over_mref,
                                        CMatrix& Stilde) {
  RealSymmetricMatrix KtildeOrInverse;
  getKtildeOrInverseFromElab(Elab_over_mref, KtildeOrInverse);
  calcCayleyTransformMatrix(KtildeOrInverse, Stilde, m_Kmat);
}

void BoxQuantization::getStildeFromEcm(double Ecm_over_mref, CMatrix& Stilde) {
  RealSymmetricMatrix KtildeOrInverse;
  getKtildeOrInverseFromEcm(Ecm_over_mref, KtildeOrInverse);
  calcCayleyTransformMatrix(KtildeOrInverse, Stilde, m_Kmat);
}

//  computes the quantization matrix and its eigenvaleus
//  specified by the "qctype" enum

void BoxQuantization::getQCMatrixFromElab(
    double Elab_over_mref, CMatrix& QC,
    QuantCondType qctype) {
  get_qc_matrix(Elab_over_mref, true, qctype, QC);
}

void BoxQuantization::getQCMatrixFromEcm(
    double Ecm_over_mref, CMatrix& QC,
    QuantCondType qctype) {
  get_qc_matrix(Ecm_over_mref, false, qctype, QC);
}

std::vector<cmplx>
BoxQuantization::getQCEigenvaluesFromElab(double Elab_over_mref,
                                          QuantCondType qctype) {
  return get_qc_eigenvalues(Elab_over_mref, true, qctype);
}

std::vector<cmplx>
BoxQuantization::getQCEigenvaluesFromEcm(double Ecm_over_mref,
                                         QuantCondType qctype) {
  return get_qc_eigenvalues(Ecm_over_mref, false, qctype);
}

//  computes Omega(mu, Q) where Q is the quantization
//  matrix specified by the "qctype" enum

cmplx BoxQuantization::getOmegaFromElab(double mu, double Elab_over_mref,
                                        QuantCondType qctype) {
  return get_omega(mu, Elab_over_mref, true, qctype);
}

cmplx BoxQuantization::getOmegaFromElab(double mu, double Elab_over_mref,
                                       const ComplexHermitianMatrix& B,
                                       QuantCondType qctype) {
  return get_omega(mu, Elab_over_mref, B, true, qctype);
}

cmplx BoxQuantization::getOmegaFromEcm(double mu, double Ecm_over_mref,
                                       QuantCondType qctype) {
  return get_omega(mu, Ecm_over_mref, false, qctype);
}

cmplx BoxQuantization::getOmegaFromEcm(double mu, double Ecm_over_mref,
                                       const ComplexHermitianMatrix& B,
                                       QuantCondType qctype) {
  return get_omega(mu, Ecm_over_mref, B, false, qctype);
}

std::list<double>
BoxQuantization::getFreeTwoParticleEnergiesInElab(double min_Elab_over_mref,
                                                  double max_Elab_over_mref) const
{
  std::list<double> all_energies;                  // final container

  auto box_it = m_boxes.begin();
  for (std::size_t chan = 0; chan < m_decay_channel_S_num.size(); ++chan) {
    // Collect NI energies for the *first* BoxMatrix in this channel
    std::list<double> tmp =
        box_it->first->getEcmTransform()
                     .getAllFreeTwoParticleEnergiesInElab(min_Elab_over_mref,
                                                       max_Elab_over_mref);
    all_energies.splice(all_energies.end(), std::move(tmp));   // O(1)

    // Skip the remaining S–blocks that belong to this channel
    std::advance(box_it, m_decay_channel_S_num[chan]);
  }

  // sort and remove duplicates
  all_energies.sort();
  all_energies.unique();

  return all_energies;
}

list<double>
BoxQuantization::getFreeTwoParticleEnergiesInEcm(double min_Ecm_over_mref,
                                            double max_Ecm_over_mref) const {
  std::list<double> all_energies;                  // final container

  auto box_it = m_boxes.begin();
  for (std::size_t chan = 0; chan < m_decay_channel_S_num.size(); ++chan) {
    // Collect NI energies for the *first* BoxMatrix in this channel
    std::list<double> tmp =
        box_it->first->getEcmTransform()
                     .getAllFreeTwoParticleEnergiesInEcm(min_Ecm_over_mref,
                                                       max_Ecm_over_mref);
    all_energies.splice(all_energies.end(), std::move(tmp));   // O(1)

    // Skip the remaining S–blocks that belong to this channel
    std::advance(box_it, m_decay_channel_S_num[chan]);
  }

  // sort and remove duplicates
  all_energies.sort();
  all_energies.unique();

  return all_energies;
}

void BoxQuantization::getEcmRootsInEcmInterval(double mu, double Ecm_over_mref_min,
                                            double Ecm_over_mref_max,
                                            QuantCondType qctype,
                                            AdaptiveBracketConfig P,
                                            std::vector<double>& roots,
                                            std::vector<uint>& fn_calls) {
  get_roots_in_interval_bracketed_by_NIs(
    mu, Ecm_over_mref_min, Ecm_over_mref_max, false, qctype, P, roots,
    fn_calls);
}


// fn calls is the number of function calls made for each NI/endpoint interval
void BoxQuantization::getEcmRootsInElabInterval(
    double mu, double Elab_over_mref_min, double Elab_over_mref_max,
    QuantCondType qctype, AdaptiveBracketConfig P, std::vector<double>& roots,
    std::vector<uint>& fn_calls) {
  get_roots_in_interval_bracketed_by_NIs(
    mu, Elab_over_mref_min, Elab_over_mref_max, true, qctype, P, roots,
    fn_calls);
}

// -----------------------------------------------------------------------------
// predict the energy shifts in a user-supplied energy window, using the shift
// observations and non-interacting pairs (NIs) as guidance
// -----------------------------------------------------------------------------
void BoxQuantization::getDeltaEcmPredictionsInEcmInterval(
        double  mu,
        double  Ecm_over_mref_min,
        double  Ecm_over_mref_max,
        QuantCondType          qctype,
        const AdaptiveBracketConfig  P,
        const std::vector<std::pair<double, NonInteractingPair>>& shift_obs_w_NIs,
        std::vector<double>&   shift_predictions,
        std::vector<uint>&     fn_calls,
        double                 guard_tol_frac) {
  get_DeltaE_predictions(
      mu,
      Ecm_over_mref_min,
      Ecm_over_mref_max,
      /* Elab = */ false,
      qctype,
      P,
      shift_obs_w_NIs,
      shift_predictions,
      fn_calls,
      guard_tol_frac);
}

void BoxQuantization::getDeltaElabPredictionsInElabInterval(
        double  mu,
        double  Elab_over_mref_min,
        double  Elab_over_mref_max,
        QuantCondType          qctype,
        const AdaptiveBracketConfig  P,
        const std::vector<std::pair<double, NonInteractingPair>>& shift_obs_w_NIs,
        std::vector<double>&   shift_predictions,
        std::vector<uint>&     fn_calls,
        double                 guard_tol_frac) {
  get_DeltaE_predictions(
      mu,
      Elab_over_mref_min,
      Elab_over_mref_max,
      /* Elab = */ true,
      qctype,
      P,
      shift_obs_w_NIs,
      shift_predictions,
      fn_calls,
      guard_tol_frac);
}

/**
 * @brief Optimized energy shift prediction working directly in Ecm coordinates
 * 
 * This function performs root finding directly in the center-of-mass frame to avoid 
 * redundant coordinate transformations. The quantization condition Ω(Ecm) = 0 is 
 * solved in CM frame, then energy shifts are computed in the requested output frame.
 * 
 * @param mu Filter parameter for Omega function
 * @param Ecm_over_mref_min Minimum CM energy for root finding
 * @param Ecm_over_mref_max Maximum CM energy for root finding
 * @param qctype Quantization condition type (StildeCB, KtildeB, etc.)
 * @param P Adaptive bracket configuration for root finding
 * @param shift_obs_w_NIs Vector of observed energy shifts with non-interacting pairs
 * @param shift_predictions Output vector of predicted energy shifts
 * @param fn_calls Output vector of function evaluation counts per interval
 * @param output_in_lab_frame If true, output shifts in lab frame; if false, in CM frame
 * @param guard_tol_frac Tolerance fraction to check other side of NI interval
 * 
 * @performance Optimizations implemented:
 * - Root finding performed directly in CM frame to avoid coordinate transformations
 */
void BoxQuantization::getDeltaEnergyPredictionsOptimized(
        double  mu,
        double  Ecm_over_mref_min,
        double  Ecm_over_mref_max,
        QuantCondType          qctype,
        const AdaptiveBracketConfig  P,
        const std::vector<std::pair<double, NonInteractingPair>>& shift_obs_w_NIs,
        std::vector<double>&   shift_predictions,
        std::vector<uint>&     fn_calls,
        bool                   output_in_lab_frame,
        double                 guard_tol_frac) {
  
  // ---- 1. Get non-interacting energies in CM frame -------------------------
  std::list<double> NI_energies_cm = getFreeTwoParticleEnergiesInEcm(
      Ecm_over_mref_min, Ecm_over_mref_max);
  
  if (NI_energies_cm.empty())
    throw std::runtime_error("No NI energies in requested CM interval.");
  
  std::vector<double> NI_vec_cm(NI_energies_cm.begin(), NI_energies_cm.end());
  
  // ---- 2. Process observations and convert to CM frame ------------------
  std::vector<double> Ecm_obs;                    // observed energies in CM frame
  std::vector<double> Ecm_free_obs;               // matching NI free energies in CM frame
  std::vector<NonInteractingPair> ni_pairs;       // store NI pairs for deferred lab frame calculations
  std::unordered_map<uint, std::reference_wrapper<const EcmTransform>> ecm_transforms;
  std::vector<uint> NI_interval_ids;              // intervals to search
  
  for (const auto& obs : shift_obs_w_NIs) {
    double dE_obs = obs.first;
    const NonInteractingPair& ni = obs.second;
    
    // Cache EcmTransform for this decay channel
    auto [it, inserted] = ecm_transforms.try_emplace(
        ni.decay_channel_idx,
        std::cref(getDecayChannelEcmTransform(ni.decay_channel_idx)));
    const EcmTransform& et = it->second.get();
    
    // Get free two-particle energies - only CM frame for now
    double Ecm_free = et.getFreeTwoParticleEnergyInEcm(ni.d1_sqr, ni.d2_sqr);
    
    // Convert observed energy to CM frame
    // Assumption: dE_obs is in lab frame initially. TODO: generalize
    // We need lab frame free energy for this conversion
    double Elab_free = et.getFreeTwoParticleEnergyInElab(ni.d1_sqr, ni.d2_sqr);
    double Elab_obs = dE_obs + Elab_free;
    double Ecm_obs_val = getEcmOverMrefFromElab(Elab_obs);
    
    // Check if CM energy is within bounds
    if (Ecm_obs_val < Ecm_over_mref_min || Ecm_obs_val > Ecm_over_mref_max)
      throw std::invalid_argument("Observation outside CM energy interval.");
    
    Ecm_obs.push_back(Ecm_obs_val);
    Ecm_free_obs.push_back(Ecm_free);
    ni_pairs.push_back(ni);  // Store for deferred lab frame calculations
    
    // Locate NI interval in CM frame
    uint interval = 0;
    double left = Ecm_over_mref_min;
    for (double pole : NI_vec_cm) {
      if (Ecm_obs_val <= pole) break;
      left = pole;
      ++interval;
    }
    NI_interval_ids.push_back(interval);
    
    // Apply guard tolerance for roots that cross the interval boundaries
    double right = (interval < NI_vec_cm.size()) ? NI_vec_cm[interval] : Ecm_over_mref_max;
    double width = right - left;
    double distL = Ecm_obs_val - left;
    double distR = right - Ecm_obs_val;
    
    if (width > 0) {
      if (distL / width < guard_tol_frac && interval > 0)
        NI_interval_ids.push_back(interval - 1);
      if (distR / width < guard_tol_frac && interval < NI_vec_cm.size())
        NI_interval_ids.push_back(interval + 1);
    }
  }
  
  // ---- 3. Deduplicate interval list -------------------------------------
  std::sort(NI_interval_ids.begin(), NI_interval_ids.end());
  NI_interval_ids.erase(std::unique(NI_interval_ids.begin(), NI_interval_ids.end()), 
                        NI_interval_ids.end());
  
  // ---- 4. Root finding in CM frame ----------------------------------
  std::vector<double> Ecm_pred;          // predicted root energies in CM frame
  fn_calls.clear();
  const double eps = 1e-9;
  
  for (uint interval : NI_interval_ids) {
    double Ecm_min = (interval == 0) ? Ecm_over_mref_min : NI_vec_cm[interval-1];
    double Ecm_max = (interval < NI_vec_cm.size()) ? NI_vec_cm[interval] : Ecm_over_mref_max;
    
    std::vector<double> roots;
    uint ncall = get_roots_in_interval(mu,
                                       Ecm_min + eps,
                                       Ecm_max - eps,
                                       /* Elab = */ false,  // Work in CM frame
                                       qctype,
                                       P,
                                       roots);
    fn_calls.push_back(ncall);
    Ecm_pred.insert(Ecm_pred.end(), roots.begin(), roots.end());
  }
  
  shift_predictions.resize(Ecm_obs.size(), 0.0);
  
  if (Ecm_pred.empty()) {
    std::cout << "\033[31m\u25A0\033[0m No roots found in CM frame" << std::endl;
    return;
  }
  
  // ---- 5. Match roots to observations and compute energy shifts -------
  for (std::size_t k = 0; k < Ecm_obs.size(); ++k) {
    double best_diff = std::numeric_limits<double>::max();
    std::size_t best_j = static_cast<std::size_t>(-1);
    
    // Find closest predicted root to observed energy (both in CM frame)
    for (std::size_t j = 0; j < Ecm_pred.size(); ++j) {
      double diff = std::abs(Ecm_pred[j] - Ecm_obs[k]);
      if (diff < best_diff) {
        best_diff = diff;
        best_j = j;
      }
    }
    
    if (best_j == static_cast<std::size_t>(-1))
      throw std::runtime_error("No root energy close to observation.");
    
    // Compute energy shift in requested output frame
    if (output_in_lab_frame) {
      // Convert predicted CM root to lab frame
      double Elab_pred = getElabOverMrefFromEcm(Ecm_pred[best_j]);
      
      // Get lab frame free energy for this observation (computed on demand)
      const NonInteractingPair& ni = ni_pairs[k];
      const EcmTransform& et = ecm_transforms.at(ni.decay_channel_idx);
      double Elab_free = et.getFreeTwoParticleEnergyInElab(ni.d1_sqr, ni.d2_sqr);
      
      // Energy shift in lab frame = predicted lab energy - free lab energy
      shift_predictions[k] = Elab_pred - Elab_free;
    } else {
      // Energy shift in CM frame = predicted CM energy - free CM energy
      shift_predictions[k] = Ecm_pred[best_j] - Ecm_free_obs[k];
    }
  }
}


void BoxQuantization::getDeltaERootsInEcmInterval(double mu, double Ecm_over_mref_min,
                             double Ecm_over_mref_max, QuantCondType qctype,
                             AdaptiveBracketConfig P,
                             std::vector<double>& roots,
                             std::vector<uint>& fn_calls) {
  get_deltaE_roots_in_interval_bracketed_by_NIs(
    mu, Ecm_over_mref_min, Ecm_over_mref_max, false, qctype, P, roots,
    fn_calls);
}

void BoxQuantization::getDeltaERootsInElabInterval(double mu, double Ecm_over_mref_min,
                             double Ecm_over_mref_max, QuantCondType qctype,
                             AdaptiveBracketConfig P,
                             std::vector<double>& roots,
                             std::vector<uint>& fn_calls) {
  get_deltaE_roots_in_interval_bracketed_by_NIs(
    mu, Ecm_over_mref_min, Ecm_over_mref_max, true, qctype, P, roots,
    fn_calls);
}


// assigns a complex value to an element in a complex Hermitian matrix
// which is either of type "ComplexHermitianMatrix" or "CMatrix"
// if "herm" is true, assigns to "CMh"; else to "CM"

void BoxQuantization::assign(const cmplx& value, uint row, uint col, bool herm,
                             ComplexHermitianMatrix& CMh, CMatrix& CM) {
  if (herm)
    CMh.put(row, col, value);
  else {
    CM.put(row, col, value);
    CM.put(col, row, conj(value));
  }
}

// assigns a real value to an element in a real symmetric matrix
// which is either of type "RealSymmetricMatrix" or "RMatrix"
// if "herm" is true, assigns to "RMh"; else to "RM"

void BoxQuantization::assign(double value, uint row, uint col, bool herm,
                             RealSymmetricMatrix& RMh, RMatrix& RM) {
  if (herm)
    RMh(row, col) = value;
  else {
    RM(row, col) = value;
    RM(col, row) = value;
  }
}

// main routine which computes the Hermitian Box matrix
// if "Elab" is true, computes using lab frame energies;
// if "Elab" is false, computes using cm frame energies;
// if "herm" is true, returns result in "Bh", else returns in "B"

void BoxQuantization::get_box_matrix(double E_over_mref,
                                     ComplexHermitianMatrix& Bh, CMatrix& B,
                                     bool Elab, bool herm) {
  if (Elab)
    for (list<pair<BoxMatrix*, uint>>::const_iterator it = m_boxes.begin();
         it != m_boxes.end(); it++) {
      it->first->setElementsFromElab(E_over_mref);
    }
  else
    for (list<pair<BoxMatrix*, uint>>::const_iterator it = m_boxes.begin();
         it != m_boxes.end(); it++) {
      it->first->setElementsFromEcm(E_over_mref);
    }
  uint bsize = m_basis.size();
  if (herm)
    Bh.resize(bsize);
  else
    B.resize(bsize, bsize);
  uint row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    BoxMatrix* bptr = rt->getBoxMatrixPtr();
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      if ((rt->getChannelIndex() != ct->getChannelIndex()) ||
          (rt->getStimestwo() != ct->getStimestwo()))
        assign(complex<double>(0, 0), row, col, herm, Bh, B);
      else {
        BoxMatrixQuantumNumbers bqn(rt->getJtimestwo(), rt->getL(),
                                    rt->getOccurrence(), ct->getJtimestwo(),
                                    ct->getL(), ct->getOccurrence());
        assign(bptr->getElement(bqn), row, col, herm, Bh, B);
      }
    }
  }
}

// main routine which computes the scattering Ktilde matrix or its inverse
//   if "Elab" is true, computes using lab frame energies;
//   if "Elab" is false, computes using cm frame energies;
//   if "herm" is true, returns result in "Kh", else returns in "K"
// The pointer "evalptr" determines if Ktilde or its inverse are evaluated

template <typename T>
void BoxQuantization::get_ktilde_matrix(double E_over_mref,
                                        RealSymmetricMatrix& Kh, RMatrix& K,
                                        bool Elab, bool herm, T* evalptr) {
  double Ecm = (Elab)
                   ? m_boxes.front().first->getEcmOverMrefFromElab(E_over_mref)
                   : E_over_mref;
  if (Ecm < m_boxes.front().first->getEcmTransform().getMinEcmOverMref())
    throw(std::invalid_argument("Ecm less than minimum not allowed"));
  uint bsize = m_basis.size();
  if (herm)
    Kh.resize(bsize);
  else
    K.resize(bsize, bsize);
  uint row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      if ((rt->getJtimestwo() != ct->getJtimestwo()) ||
          (rt->getOccurrence() != ct->getOccurrence()))
        assign(0.0, row, col, herm, Kh, K);
      else {
        double kres = evalptr->calculate(
            rt->getJtimestwo(), rt->getL(), rt->getStimestwo(),
            rt->getChannelIndex(), ct->getL(), ct->getStimestwo(),
            ct->getChannelIndex(), Ecm);
        assign(kres, row, col, herm, Kh, K);
      }
    }
  }
}

void BoxQuantization::assign_matrices(double E_over_mref, bool Elab,
                                      ComplexHermitianMatrix& B,
                                      RealSymmetricMatrix& KvOrInv) {
  if (Elab) {
    getBoxMatrixFromElab(E_over_mref, B);
    if (m_Kmat != 0)
      getKtildeFromElab(E_over_mref, KvOrInv);
    else
      getKtildeinvFromElab(E_over_mref, KvOrInv);
  } else {
    getBoxMatrixFromEcm(E_over_mref, B);
    if (m_Kmat != 0)
      getKtildeFromEcm(E_over_mref, KvOrInv);
    else
      getKtildeinvFromEcm(E_over_mref, KvOrInv);
  }
}

void BoxQuantization::assign_ct_matrices(double E_over_mref, bool Elab,
                                         CMatrix& CB, CMatrix& Stilde) {
  if (Elab) {
    getCTBoxMatrixFromElab(E_over_mref, CB);
    getStildeFromElab(E_over_mref, Stilde);
  } else {
    getCTBoxMatrixFromEcm(E_over_mref, CB);
    getStildeFromEcm(E_over_mref, Stilde);
  }
}

template <typename T>
set<BoxQuantBasisState>
BoxQuantization::find_excluded_states_from_ktilde(T* evalptr) {
  set<BoxQuantBasisState> exclusions;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++) {
    bool exclude = true;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++) {
      if ((rt->getJtimestwo() == ct->getJtimestwo()) &&
          (rt->getOccurrence() == ct->getOccurrence()) &&
          (!(evalptr->isZero(rt->getJtimestwo(), rt->getL(), rt->getStimestwo(),
                             rt->getChannelIndex(), ct->getL(),
                             ct->getStimestwo(), ct->getChannelIndex())))) {
        exclude = false;
        break;
      }
    }
    if (exclude)
      exclusions.insert(*rt);
  }
  return exclusions;
}

void BoxQuantization::output_matrices(double E_over_mref, bool Elab,
                                      std::ostream& fout) {
  fout.precision(15);
  ComplexHermitianMatrix B;
  RealSymmetricMatrix Kv;
  assign_matrices(E_over_mref, Elab, B, Kv);
  string Kname;
  if (Elab) {
    fout << "Elab:=" << E_over_mref << ":" << endl;
    Kname = (m_Kmat != 0) ? "Ktilde" : "Ktildeinv";
  } else {
    fout << "Ecm:=" << E_over_mref << ":" << endl;
    Kname = (m_Kmat != 0) ? "Ktilde" : "Ktildeinv";
  }
  uint row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    uint ap = rt->getChannelIndex();
    uint Spx2 = rt->getStimestwo();
    uint Jpx2 = rt->getJtimestwo();
    uint Lp = rt->getL();
    uint np = rt->getOccurrence();
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      uint a = ct->getChannelIndex();
      uint Sx2 = ct->getStimestwo();
      uint Jx2 = ct->getJtimestwo();
      uint L = ct->getL();
      uint n = ct->getOccurrence();
      fout << "B[[" << ap << "," << Spx2 << "," << Jpx2 << "," << Lp << ","
           << np << "],[" << a << "," << Sx2 << "," << Jx2 << "," << L << ","
           << n << "]]:=cmplx" << B(row, col) << ":" << endl;
    }
  }
  row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    uint ap = rt->getChannelIndex();
    uint Spx2 = rt->getStimestwo();
    uint Jpx2 = rt->getJtimestwo();
    uint Lp = rt->getL();
    uint np = rt->getOccurrence();
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      uint a = ct->getChannelIndex();
      uint Sx2 = ct->getStimestwo();
      uint Jx2 = ct->getJtimestwo();
      uint L = ct->getL();
      uint n = ct->getOccurrence();
      fout << Kname << "[[" << ap << "," << Spx2 << "," << Jpx2 << "," << Lp
           << "," << np << "],[" << a << "," << Sx2 << "," << Jx2 << "," << L
           << "," << n << "]]:=" << Kv(row, col) << ":" << endl;
    }
  }
}

void BoxQuantization::output_ct_matrices(double E_over_mref, bool Elab,
                                         std::ostream& fout) {
  fout.precision(15);
  CMatrix CB, Stilde;
  assign_ct_matrices(E_over_mref, Elab, CB, Stilde);
  uint row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    uint ap = rt->getChannelIndex();
    uint Spx2 = rt->getStimestwo();
    uint Jpx2 = rt->getJtimestwo();
    uint Lp = rt->getL();
    uint np = rt->getOccurrence();
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      uint a = ct->getChannelIndex();
      uint Sx2 = ct->getStimestwo();
      uint Jx2 = ct->getJtimestwo();
      uint L = ct->getL();
      uint n = ct->getOccurrence();
      fout << "CB[[" << ap << "," << Spx2 << "," << Jpx2 << "," << Lp << ","
           << np << "],[" << a << "," << Sx2 << "," << Jx2 << "," << L << ","
           << n << "]]:=cmplx" << CB(row, col) << ":" << endl;
    }
  }
  row = 0;
  for (set<BoxQuantBasisState>::iterator rt = m_basis.begin();
       rt != m_basis.end(); rt++, row++) {
    uint ap = rt->getChannelIndex();
    uint Spx2 = rt->getStimestwo();
    uint Jpx2 = rt->getJtimestwo();
    uint Lp = rt->getL();
    uint np = rt->getOccurrence();
    uint col = 0;
    for (set<BoxQuantBasisState>::iterator ct = m_basis.begin();
         ct != m_basis.end(); ct++, col++) {
      uint a = ct->getChannelIndex();
      uint Sx2 = ct->getStimestwo();
      uint Jx2 = ct->getJtimestwo();
      uint L = ct->getL();
      uint n = ct->getOccurrence();
      fout << "Stilde[[" << ap << "," << Spx2 << "," << Jpx2 << "," << Lp << ","
           << np << "],[" << a << "," << Sx2 << "," << Jx2 << "," << L << ","
           << n << "]]:=" << Stilde(row, col) << ":" << endl;
    }
  }
}

void BoxQuantization::get_qc_matrix(double E_over_mref, bool Elab,
                                    QuantCondType qctype, CMatrix& Q) {
  if (qctype == QuantCondType::StildeCB) {
    CMatrix CB, Stilde;
    assign_ct_matrices(E_over_mref, Elab, CB, Stilde);
    int N = CB.size(0);
    Q.resize(N, N);
    cmplx one(1.0, 0.0);
    // evaluate Q = 1 + Stilde * CB
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        cmplx z(0.0, 0.0);
        for (int m = 0; m < N; ++m) {
          z += Stilde(row, m) * CB(m, col);
        }
        if (row == col)
          z += one;
        Q.put(row, col, z);
      }
  } else if (qctype == QuantCondType::StildeinvCB) {
    CMatrix CB, Stilde;
    assign_ct_matrices(E_over_mref, Elab, CB, Stilde);
    int N = CB.size(0);
    Q.resize(N, N);
    // evaluate Q = Stildeinv + CB
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        Q.put(row, col, std::conj(Stilde(col, row)) + CB(row, col));
      }
  } else if (qctype == QuantCondType::KtildeB) {
    if (isKtildeInverseMode()) {
      throw(std::invalid_argument("Cannot compute 1-Ktilde*B quantization "
                                  "condition in KtildeInverse mode"));
    }
    ComplexHermitianMatrix B;
    RealSymmetricMatrix Ktilde;
    assign_matrices(E_over_mref, Elab, B, Ktilde);
    int N = B.size();
    Q.resize(N, N);
    cmplx one(1.0, 0.0);
    // evaluate Q = 1 - Ktilde * B
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        cmplx z(0.0, 0.0);
        for (int m = 0; m < N; ++m) {
          z -= Ktilde(row, m) * B(m, col);
        }
        if (row == col)
          z += one;
        Q.put(row, col, z);
      }
  } else if (qctype == QuantCondType::KtildeinvB) {
    if (isKtildeMode()) {
      throw(std::invalid_argument(
          "Cannot compute Ktildeinv-B quantization condition in Ktilde mode"));
    }
    ComplexHermitianMatrix B;
    RealSymmetricMatrix Ktildeinv;
    assign_matrices(E_over_mref, Elab, B, Ktildeinv);
    int N = B.size();
    Q.resize(N, N);
    // evaluate Q = Ktildeinv - B
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        Q.put(row, col, cmplx(Ktildeinv(row, col), 0.0) - B(row, col));
      }
  }
}

void BoxQuantization::get_qc_matrix(double E_over_mref, 
                                    const ComplexHermitianMatrix& B,
                                    bool Elab, QuantCondType qctype,
                                    CMatrix& Q) {
  if (qctype == QuantCondType::StildeCB) {
    // Compute Cayley transform of B: CB = (1 + i*B) * (1 - i*B)^(-1)
    CMatrix CB;
    calcCayleyTransformMatrix(B, CB, true);
    // Compute Stilde from K-matrix
    CMatrix Stilde;
    RealSymmetricMatrix KtildeOrInverse;
    if (Elab) {
      getKtildeOrInverseFromElab(E_over_mref, KtildeOrInverse);
    } else {
      getKtildeOrInverseFromEcm(E_over_mref, KtildeOrInverse);
    }
    calcCayleyTransformMatrix(KtildeOrInverse, Stilde, m_Kmat);
    
    int N = CB.size(0);
    Q.resize(N, N);
    cmplx one(1.0, 0.0);
    // evaluate Q = 1 + Stilde * CB
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        cmplx z(0.0, 0.0);
        for (int m = 0; m < N; ++m) {
          z += Stilde(row, m) * CB(m, col);
        }
        if (row == col)
          z += one;
        Q.put(row, col, z);
      }
  } else if (qctype == QuantCondType::StildeinvCB) {
    // Compute Cayley transform of B: CB = (1 + i*B) * (1 - i*B)^(-1)
    CMatrix CB;
    calcCayleyTransformMatrix(B, CB, true);
    // Compute Stilde from K-matrix
    CMatrix Stilde;
    RealSymmetricMatrix KtildeOrInverse;
    if (Elab) {
      getKtildeOrInverseFromElab(E_over_mref, KtildeOrInverse);
    } else {
      getKtildeOrInverseFromEcm(E_over_mref, KtildeOrInverse);
    }
    calcCayleyTransformMatrix(KtildeOrInverse, Stilde, m_Kmat);
    
    int N = CB.size(0);
    Q.resize(N, N);
    // evaluate Q = Stildeinv + CB
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        Q.put(row, col, std::conj(Stilde(col, row)) + CB(row, col));
      }
  } else if (qctype == QuantCondType::KtildeB) {
    if (isKtildeInverseMode()) {
      throw(std::invalid_argument("Cannot compute 1-Ktilde*B quantization "
                                  "condition in KtildeInverse mode"));
    }
    // Get Ktilde matrix
    RealSymmetricMatrix Ktilde;
    if (Elab) {
      getKtildeFromElab(E_over_mref, Ktilde);
    } else {
      getKtildeFromEcm(E_over_mref, Ktilde);
    }
    
    int N = B.size();
    Q.resize(N, N);
    cmplx one(1.0, 0.0);
    // evaluate Q = 1 - Ktilde * B
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        cmplx z(0.0, 0.0);
        for (int m = 0; m < N; ++m) {
          z -= Ktilde(row, m) * B(m, col);
        }
        if (row == col)
          z += one;
        Q.put(row, col, z);
      }
  } else if (qctype == QuantCondType::KtildeinvB) {
    if (isKtildeMode()) {
      throw(std::invalid_argument(
          "Cannot compute Ktildeinv-B quantization condition in Ktilde mode"));
    }
    // Get Ktildeinv matrix
    RealSymmetricMatrix Ktildeinv;
    if (Elab) {
      getKtildeinvFromElab(E_over_mref, Ktildeinv);
    } else {
      getKtildeinvFromEcm(E_over_mref, Ktildeinv);
    }
    
    int N = B.size();
    Q.resize(N, N);
    // evaluate Q = Ktildeinv - B
    for (int row = 0; row < N; ++row)
      for (int col = 0; col < N; ++col) {
        Q.put(row, col, cmplx(Ktildeinv(row, col), 0.0) - B(row, col));
      }
  }
}

std::vector<cmplx> BoxQuantization::get_qc_eigenvalues(double E_over_mref,
                                                       bool Elab,
                                                       QuantCondType qctype) {
  CMatrix Q;
  get_qc_matrix(E_over_mref, Elab, qctype, Q);
  Cvector eigvals;
  Diagonalizer D;
  D.getEigenvalues(Q, eigvals);
  return eigvals;
}

cmplx BoxQuantization::get_omega(double mu, double E_over_mref, bool Elab,
                                 QuantCondType qctype) {
  CMatrix Q;
  get_qc_matrix(E_over_mref, Elab, qctype, Q);
  DeterminantCalculator DC;
  if (mu >= 0.0) {
    return DC.getOmega(mu, Q);
  }
  // mu < 0.0, dont calculate omega
  return DC.getDeterminant(Q);
}

// Version that takes an explicit B matrix instead of computing it internally
cmplx BoxQuantization::get_omega(double mu, double E_over_mref,
                                  const ComplexHermitianMatrix& B,
                                  bool Elab,
                                  QuantCondType qctype) {
  CMatrix Q;
  get_qc_matrix(E_over_mref, B, Elab, qctype, Q);
  DeterminantCalculator DC;
  if (mu >= 0.0) {
    return DC.getOmega(mu, Q);
  }
  // mu < 0.0, dont calculate omega
  return DC.getDeterminant(Q);
}

// -----------------------------------------------------------------------------
//  Predict shifts by matching on raw E_lab/cm values first
// -----------------------------------------------------------------------------
void BoxQuantization::get_DeltaE_predictions(
        double                      mu,
        double                      E_over_mref_min,
        double                      E_over_mref_max,
        bool                        Elab,
        QuantCondType               qctype,
        const AdaptiveBracketConfig       P,
        const std::vector<std::pair<double, NonInteractingPair>>& shift_obs_w_NIs,
        std::vector<double>&        shift_predictions,   // output
        std::vector<uint>&          fn_calls,            // omega evals per interval
        double                      guard_tol_frac) {
  // ---- 1.  NI intervals -------------------------------------------------------
  std::list<double> NI_energies = Elab
        ? getFreeTwoParticleEnergiesInElab(E_over_mref_min, E_over_mref_max)
        : getFreeTwoParticleEnergiesInEcm (E_over_mref_min, E_over_mref_max);
  if (NI_energies.empty())
      throw std::runtime_error("No NI energies in requested interval.");
  std::vector<double> NI_vec(NI_energies.begin(), NI_energies.end());

  // ---- 2.  Observation bookkeeping ---------------------------------------
  std::vector<double> E_obs;                     // raw energies
  std::vector<double> E_free_obs;                // matching NI free energies
  std::unordered_map<uint, std::reference_wrapper<const EcmTransform>> ecm_transforms;    // per channel cache
  std::vector<uint> NI_interval_ids;             // intervals to search

  for (const auto& obs : shift_obs_w_NIs)
  {
    double dE_obs = obs.first;
    const NonInteractingPair& ni = obs.second;

    // EcmTransform cache
    auto [it, inserted] = ecm_transforms.try_emplace(
      ni.decay_channel_idx,
      std::cref(getDecayChannelEcmTransform(ni.decay_channel_idx)));
    const EcmTransform& et = it->second.get();

    // free two-particle energy for this datum
    double E_free = Elab
          ? et.getFreeTwoParticleEnergyInElab(ni.d1_sqr, ni.d2_sqr)
          : et.getFreeTwoParticleEnergyInEcm(ni.d1_sqr, ni.d2_sqr);

    double E = dE_obs + E_free;
    if (E < E_over_mref_min || E > E_over_mref_max)
        throw std::invalid_argument("Shift observation outside interval.");

    E_obs.push_back(E);
    E_free_obs.push_back(E_free);

    // locate NI interval
    uint interval = 0;
    double left   = E_over_mref_min;
    for (double pole : NI_vec) {
        if (E <= pole) break;
        left = pole;
        ++interval;
    }
    NI_interval_ids.push_back(interval);

    // guard if near a pole
    double right  = (interval < NI_vec.size()) ? NI_vec[interval]
                                               : E_over_mref_max;
    double width  = right - left;
    double distL  = E - left;
    double distR  = right - E;
    if (width > 0) {
        if (distL / width < guard_tol_frac && interval > 0)
            NI_interval_ids.push_back(interval - 1);
        if (distR / width < guard_tol_frac &&
            interval < NI_vec.size())
            NI_interval_ids.push_back(interval + 1);
    }
  }

  // ---- 3.  Deduplicate interval list -------------------------------------
  std::sort(NI_interval_ids.begin(), NI_interval_ids.end());
  NI_interval_ids.erase(std::unique(NI_interval_ids.begin(),
                                    NI_interval_ids.end()),
                        NI_interval_ids.end());

  // ---- 4.  Root finding ---------------------------------------------------
  std::vector<double> E_pred;          // root energies (Ecm or Elab)
  fn_calls.clear();
  const double eps = 1e-9;

  for (uint interval : NI_interval_ids) {
    double Emin = (interval == 0) ? E_over_mref_min : NI_vec[interval-1];
    double Emax = (interval < NI_vec.size()) ? NI_vec[interval]
                                             : E_over_mref_max;

    std::vector<double> roots;
    uint ncall = get_roots_in_interval(mu,
                                       Emin + eps,
                                       Emax - eps,
                                       Elab,
                                       qctype,
                                       P,
                                       roots);
    fn_calls.push_back(ncall);
    E_pred.insert(E_pred.end(), roots.begin(), roots.end());
  }
  shift_predictions.resize(E_obs.size(), 0.0);

  if (E_pred.empty()) {
    cout << "\033[31m\u25A0\033[0m No roots found" << endl;
    return;
  }

  // // ---- 5.  Greedy nearest-energy matching --------------------------
  // shift_predictions.resize(E_obs.size(),
  //                          std::numeric_limits<double>::quiet_NaN());

  for (std::size_t k = 0; k < E_obs.size(); ++k)
  {
    double best_diff = std::numeric_limits<double>::max();
    std::size_t best_j = static_cast<std::size_t>(-1);

    for (std::size_t j = 0; j < E_pred.size(); ++j) {
      double diff = std::abs(E_pred[j] - E_obs[k]);
      if (diff < best_diff) { best_diff = diff; best_j = j; }
    }

    if (best_j == static_cast<std::size_t>(-1))
      throw std::runtime_error("No root energy close to observation.");

    // Delta E_pred = matched root − datum's own free NI energy
    shift_predictions[k] = E_pred[best_j] - E_free_obs[k];
  }
}


void BoxQuantization::get_deltaE_roots_in_interval_bracketed_by_NIs(
      double mu, double E_over_mref_min, double E_over_mref_max, bool Elab,
      QuantCondType qctype, AdaptiveBracketConfig P, std::vector<double>& roots,
      std::vector<uint>& fn_calls) {
  std::vector<double> intervals;
  list<double> NIs;
  intervals.push_back(E_over_mref_min);
  if (Elab) {
    NIs = getFreeTwoParticleEnergiesInElab(E_over_mref_min, E_over_mref_max);
  }
  else {
    NIs = getFreeTwoParticleEnergiesInEcm(E_over_mref_min, E_over_mref_max);
  }
  for (list<double>::const_iterator it = NIs.begin(); it != NIs.end(); ++it) {
    intervals.push_back(*it);
  }
  intervals.push_back(E_over_mref_max);

  uint num_of_expected_roots = roots.size();

  if (NIs.size() != num_of_expected_roots) {
    throw(std::invalid_argument("Number of expected roots in Omega does not match number of NIs. "
                                "Consider adjusting the energy interval."));
  }

  std::vector<double> temp_roots;
  uint num_of_roots_found = 0;
  for (size_t i = 0; i < intervals.size() - 1; ++i) {
    temp_roots.clear();

    double E_min = intervals[i] - 1e-9;
    double E_max = intervals[i + 1] - 1e-9;
    // above might be hitting NIs and getting nan's
    // which is not good. Likely is an issue with first calculating the box
    // matrix then relying on numerical limits to work out.

    if (E_min > E_max)
      throw(std::invalid_argument("Bad interval in getRootsInEInterval"));
    // get the roots in the interval
    fn_calls.push_back(
        get_roots_in_interval(mu, E_min, E_max, Elab, qctype, P, temp_roots));
    // add the roots to the final list
    for (std::vector<double>::const_iterator it = temp_roots.begin();
         it != temp_roots.end(); ++it) {
      roots[num_of_roots_found++] = *it;
      if (num_of_roots_found == num_of_expected_roots) {
        break; // stop if we have found all expected roots
      }
    }
  }



  // now set roots to energy differences in place
  uint energy_count = 0;
  for (list<double>::const_iterator it = NIs.begin(); it != NIs.end(); ++it) {
    roots[energy_count++] -= *it; // this is the delta E
  }
}

void BoxQuantization::get_roots_in_interval_bracketed_by_NIs(double mu, double E_over_mref_min,
                                            double E_over_mref_max, bool Elab,
                                            QuantCondType qctype,
                                            AdaptiveBracketConfig P,
                                            std::vector<double>& roots,
                                            std::vector<uint>& fn_calls) {
  std::vector<double> intervals;
  list<double> NIs;
  intervals.push_back(E_over_mref_min);
  if (Elab) {
    NIs = getFreeTwoParticleEnergiesInElab(E_over_mref_min, E_over_mref_max);
  }
  else {
    NIs = getFreeTwoParticleEnergiesInEcm(E_over_mref_min, E_over_mref_max);
  }
  for (list<double>::const_iterator it = NIs.begin(); it != NIs.end(); ++it) {
    intervals.push_back(*it);
  }
  intervals.push_back(E_over_mref_max);

  for (size_t i = 0; i < intervals.size() - 1; ++i) {
    std::vector<double> temp_roots;
    double E_min = intervals[i] - 1e-9;
    double E_max = intervals[i + 1] - 1e-9;
    // above might be hitting NIs and getting nan's
    // which is not good. Likely is an issue with first calculating the box
    // matrix then relying on numerical limits to work out.

    if (E_min > E_max)
      throw(std::invalid_argument("Bad interval in getRootsInEInterval"));
    // get the roots in the interval
    fn_calls.push_back(
        get_roots_in_interval(mu, E_min, E_max, Elab, qctype, P, temp_roots));
    // add the roots to the final list
    for (std::vector<double>::const_iterator it = temp_roots.begin();
         it != temp_roots.end(); ++it) {
      if (std::find(roots.begin(), roots.end(), *it) == roots.end()) {
        roots.push_back(*it);
      }
    }
  }
}

// return the number of omega calls made for each NI/endpoint interval
uint BoxQuantization::get_roots_in_interval(double mu, double E_over_mref_min,
                                            double E_over_mref_max, bool Elab,
                                            QuantCondType qctype,
                                            AdaptiveBracketConfig P,
                                            std::vector<double>& roots) {
  std::function<cmplx(double)> f = [this, mu, Elab,
                                    qctype](double E_over_mref) {
    cmplx omega = get_omega(mu, E_over_mref, Elab, qctype);
    return omega;
  };

  AdaptiveBracketRootFinder RF(P, f);

  RF.findRoots(E_over_mref_min, E_over_mref_max, roots);
  return RF.evalCount();
}

// ***************************************************************************************
