/*!
\page Ms1Ms2MatchingExample A Mass Spectrometry Example: Matching high- and low-accuracy precursors

\section setup Setup
Another common preprocessing step in mass spectrometry data analysis is the
replacement of low-accuracy precursor masses (which come from the instrument
prescan) with their high mass accuracy couterparts (which are acquired in the
MS1 full scan). In the example at hand, we assume that the data stem from a
hybrid Thermo LTQ/Orbitrap (classic) instrument, with a prescan resolution of
7500 and a MS1 full scan resolution of 60000. 

\section method Method
We use an approach similar similar to what is sketched in the first MS example
to generate extracted ion currents (XICs). Given the XICs and all precursor
masses, we would like to determine the closest cross-set neighbors in the XIC
and precursor sets, taking into consideration the error bounds that are
associated with (i) the low-accuracy precursor, (ii) the XIC m/z, and (iii) the
XIC retention time measurements. Methodologically,

\li we use fast box intersection to derive sets of potential (XIC, parent mass)
    pairs
\li within these sets we need to select a pair by e.g. minimizing the m/z
    distance between the XIC and the precursor. In the example at hand,
    conflict resolution is kept overly simple and we accept the first
    candidate.

\section setting_up Setting Things Up
First order of business is the definition of the user classes: a class to
represent an XIC and a class to hold the scan precursor information.

\code
struct Xic
{
    double rt_;
    double mz_;
    double abundance_;
    Xic(const double &rt, const double& mz, const double & abundance) 
      : rt_(rt), mz_(mz), abundance_(abundance) {}
};

struct MS2Scan
{
    double rt_;
    double mz_;
    // some more definitions for the ions etc...
    MS2Scan(const double &rt, const double& mz) : rt_(rt), mz_(mz) {}
};
\endcode

Then we need to provide some type information to \c libfbi:

\code
namespace fbi {
template<> struct Traits<Xic> : mpl::TraitsGenerator<double, double> {};
template<> struct Traits<MS2Scan> : mpl::TraitsGenerator<double, double> {};
} 
\endcode
This tells \t libfbi that both, the XIC and the ms2scan, can be queried in two
dimensions, each of which is of type double (i.e. a floating point number).
Note, that it is not necessary to inform \c libfbi about all data members of a
class if these data members are not used in the fast box intersection process.
For example, the \c Xic struct has an \c abundance member that is irrelevant
for the intersection task at hand and is thus not included in \c Traits<Xic>.

Now we need to write the \c Xic and \c MS2Scan accessor/adapter functions that
allow libfbi to derive interval information in each dimension for each ms2scan
and XIC measurement. For the \c Xic struct, this yields

\code
struct XicBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<Xic>::key_type>::type 
  get(const Xic &) const;

  double fullScanPpm_;
  double rtWindow_;

  XicBoxGenerator(double fullScanPpm, double rtWindow)
    : fullScanPpm_(fullScanPpm), rtWindow_(rtWindow)
  {}
};

template <>
std::pair<double, double>  
XicBoxGenerator::get<0>(const Xic& xic) const 
{
  return std::make_pair(
    xic.mz_* (1 - fullScanPpm_ * 1E-6), 
    xic.mz_* (1 + fullScanPpm_ * 1E-6));
}

template <>
std::pair<double, double>  
XicBoxGenerator::get<1>(const Xic & xic) const
{
  return std::make_pair(xic.rt_ - rtWindow_, xic.rt_ + rtWindow_);
}
\endcode

And for the \c MS2Scan struct, we have

\code
struct MS2ScanBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<MS2Scan>::key_type>::type 
  get(const MS2Scan &) const;

  double preScanPpm_;
  double rtWindow_;

  MS2ScanBoxGenerator(double preScanPpm, double rtWindow)
    : preScanPpm_(preScanPpm), rtWindow_(rtWindow)
  {}
};

template <>
std::pair<double, double>  
MS2ScanBoxGenerator::get<0>(const MS2Scan& ms2scan) const 
{
  return std::make_pair(
    ms2scan.mz_* (1 - preScanPpm_ * 1E-6), 
    ms2scan.mz_* (1 + preScanPpm_ * 1E-6));
}

template <>
std::pair<double, double>  
MS2ScanBoxGenerator::get<1>(const MS2Scan & ms2scan) const
{
  return std::make_pair(ms2scan.rt_ - rtWindow_, ms2scan.rt_ + rtWindow_);
}
\endcode

We also define two helper functions, that allow us to read \c Xic and \c
MS2Scan data from text files:

\code
std::vector<Xic> parseXicFile(ProgramOptions& options);
std::vector<MS2Scan> parseMS2ScanFile(ProgramOptions & options);
\endcode

The functions that a \c ProgramOptions reference that holds program parameters;
see the implementation in the \c examples directory for their definition and
the program option parsign code (based on boost::program_options).

\section main The Main Program

With the above preparations carried out, writing the code that intersects the
XIC and precursor sets is straightforward. We start and read the data:

\code
int main(int argc, char* argv[])
{
    using namespace fbi;

    ProgramOptions options;
    if (parseProgramOptions(argc, argv, options) != 0) {
        return -1;
    }
    // load data
    std::vector<Xic> xics = parseXicFile(options);
    std::vector<MS2Scan> ms2scans = parseMS2ScanFile(options);
\endcode

With the data available, we carry out the intersection:

\code
    auto adjList = SetA<Xic, 0, 1>::SetB<MS2Scan, 0, 1>::intersect(
      xics, XicBoxGenerator(options.fullscanPpm_, options.rtWindow_),
      ms2scans, MS2ScanBoxGenerator(options.prescanPpm_, options.rtWindow_));
\endcode

The \c XicBoxGenerator and \c MS2ScanBoxGenerator instances will always return
the proper intervals for a given XIC and/or precursor. Determining the
high-resolution precursor mass for each low-resolution precursor now boils down
to resolving the conflict between the pair candidates. We choose the precursor
with the smallest m/z distance:

\code
    std::ofstream ofs(options.outputFileName_.c_str());
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    size_t nXics = xics.size();
    for (size_t i = 0; i < nXics; ++i) {
        if (!adjList[i].empty()) {
            typedef std::set<unsigned int>::const_iterator SI;
            SI best = adjList[i].begin();
            double bestDist = std::numeric_limits<double>::max();
            // resolve conflict by choosing the closest m/z
            for (SI j = adjList[i].begin(); j != adjList[i].end(); ++j) {
                size_t k = *j - nXics;
                double dist = std::abs(xics[i].mz_ - ms2scans[k].mz_);
                if (dist < bestDist) {
                    best = j;
                    bestDist = dist;
                }
            }
            // print the pair (current precursor -> new precursor)
            size_t k = *best - nXics;
            ofs << xics[i].mz_ << '\t' << xics[i].rt_
              << "\t->\t" << ms2scans[k].mz_ << '\t' << ms2scans[k].rt_ << '\n';
       }
    }
    return 0;
}
\endcode
*/
