#include <algorithm>
#include <iostream>

struct Xic {
  double mz_;
  double rt_;
};

std::ostream & operator<< ( std::ostream & os, const Xic & xic) {
  os << xic.mz_ << "\t" << xic.rt_;
  return os;
}

template <class CentroidIterType, class LabelIterType>
std::vector<Xic>
createXicVector(CentroidIterType CIT, CentroidIterType CEND, LabelIterType LIT, LabelIterType LEND, std::vector<unsigned int> & counter) {

  std::vector<Xic> xics(counter.size());
  std::fill(counter.begin(), counter.end(), 0);

  for (; CIT != CEND && LIT != LEND; ++CIT, ++LIT) {
    xics[*LIT-1].mz_ += CIT->mz_;
    xics[*LIT-1].rt_ += CIT->rt_;
    counter[*LIT-1] += 1;
  }
  for (std::vector<Xic>::size_type i = 0; i < xics.size(); ++i) {
    xics[i].mz_ /= counter[i];
    xics[i].rt_ /= counter[i];
  }
  
  return xics;
}



