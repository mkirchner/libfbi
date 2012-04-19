/*!
\page XicConstructionExample A Mass Spectrometry Example: Calculating Extracted Ion Currents (XICs)
\section setup Setup
A common problem in mass spectrometry data analysis is the preprocessing of
high-resolution liquid chromatography/mass spectrometry data (LC/MS). In
general, such data consists of a set of mass spectra that are acquired across a
liquid chromatography gradient, yielding approximately one MS1 spectrum per
second. The preprocessing task is to extract the ion current for each observed
species from these data.

In the case at hand, we follow a somewhat established procedure to extract the
XIC: 
\li We first determine the centroid of each local maximum in the m/z
domain. This step and the underlying methods have thoroughly been described in a
number of publications (e.g. [Cox and Mann, 2008]) and the remainder of the
example assumes that the data is available in the form of a list of centroids
that (to keep things simple) can be read from disk.
\li We then determine which centroids belong together to form a XIC. This is
where \c libfbi comes into play: two centroids belong to the same XIC if their
uncertainty boxes overlap. The centroid uncertainty is generally obtained by
the centroiding procedure, however, to keep things simple, we here only use an
m/z-dependent parts-per-million (ppm) estimate that allows for a simplified
description of the mass accuracy behavior of an LTQ Orbitrap instrument.
Furthermore, instead of using the retention time measurement provided by the
mass analyzer, we rely on the spectrum number; this information is more robust
as the retention time differences between runs can vary due to varying numbers
of fragementation scans that are acquired in parallel.

All of the following material is provided as the \c example-xic-construction 
sample application which is provided together with the \c libfbi package. The
complete source and header files reside in the \c example/ subdirectory of the 
\c libfbi project.

\section prep Necessary Preparations and Definitions

In order to tackle the XIC extraction task with \c libfbi, we first need to
define a data type that holds the centroid information:

\code
struct Centroid 
{
  double rt_;
  double mz_;
  double sn_;
  double abundance_;
  Centroid(const double &rt, const double& mz,
    const double & sn, const double & abundance) 
    : rt_(rt), mz_(mz), sn_(sn), abundance_(abundance) {}
};
\endcode

The \c Centrod struct holds all information available for a single centroid. We
now procedd to provide the code necessary to enable \c libfbi to derive the
uncertainty information for a centroid.

First, we tell \c libfbi that \c Centroid is a two-dimensional
measurement (we will make use of the mass/charge ratio and the spectrum number):

\code
namespace fbi {

template<>
struct Traits<Centroid> : mpl::TraitsGenerator<
  double, double> {};

} // end namespace fbi
\endcode
Note that \c Traits<Centroid> needs to be specialized in the \c fbi namespace. 

We now need to provide a wrapper/adapter class, through which \c libfbi will
access the information in the \c Centroid objects. In fact, \c libfbi never
accesses the \c Centroid objects themselves but only requests their uncertainty
ranges. In other words, we need to provide a wrapper that transforms each \c
Centroid coordinate (the m/z measurement and the spectrum number in our case)
into a \c std::pair that holds the corresponding uncertainty range.
Practically, this translates to providing an adapter class that has a template
member function \c get<Dim>() and a corresponding specialization for each
dimension that will be queried by \c libfbi. A potential adapter class looks
like this:

\code
struct BoxGenerator
{
  template <size_t N>
  typename boost::tuples::element<N, 
    typename fbi::Traits<Centroid>::key_type>::type 
  get(const Centroid &) const;

  double mzOffset_;
  double mzWindowPpm_;
  double snOffset_;
  double snWindow_;

  BoxGenerator(double mzWindowPpm, double snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm),
      snOffset_(0.0), snWindow_(snWindow)
  {}

  BoxGenerator(double mzOffset, double mzWindowPpm, 
    double snOffset, double snWindow)
    : mzOffset_(mzOffset), mzWindowPpm_(mzWindowPpm),
      snOffset_(snOffset), snWindow_(snWindow)
  {}
};
\endcode
It is recommended to use the return type specification for \c
get<T>() as illustrated in this example.

\c BoxGenerator defines a class that holds all functions and additional
information (such as the ppm values used to determine the m/z interval size, and
the allowed number of scans in which a centroid has been missed) neccessary to
derive the m/z and spectrum number intervals. The complete template
specialization of the respective interval calculation functions is thus very
simple:

\code
template <>
std::pair<double, double>  
BoxGenerator::get<0>(const Centroid & centroid) const 
{
  return std::make_pair(
    mzOffset_ + centroid.mz_* (1-  mzWindowPpm_* 1E-6), 
    mzOffset_ + centroid.mz_* (1+ mzWindowPpm_ * 1E-6) );
}

template <>
std::pair<double, double>  
BoxGenerator::get<1>(const Centroid & centroid) const
{
  return std::make_pair(
    snOffset_ + centroid.sn_ - snWindow_ - 0.3, 
    snOffset_ + centroid.sn_ + snWindow_ + 0.3 );
}
\endcode

This finalizes the setup and preparation work that needs to be carried out to
make use of \c libfbi for XIC determination. Hence, let's turn to writing the
respective main program.

\section main Using libfbi in your main program

The first step to XIC determination is to read the centroid data from disk and
into \c Centroid instances:

\code
int main(int argn, char* argv[]) {
    //
    // ... 
    //

    // read centroid data from file
    std::vector<Centroid> centroids = parseFile(options);

    //
    // ... 
    //
\endcode

XIC calculation is a self-intersection problem: within one set of boxes (each
with an (m/z, spectrum number) measurement pair at its center), we would like to
find all boxes that intersect with each other. For this kind of problem 
\c libfbi provides the \c SetA<T,Dim..>::intersect<T,F..>(...) function:

\code
    auto adjList = SetA<Centroid, 0, 1>::
      intersect(centroids, BoxGenerator(2, 2.1), BoxGenerator(2, 2.1));
\endcode

This tells \c libfbi to carry out a box intersection in dimensions 0 and 1,
using a 4ppm x 4.2 scans interval for each measurement (first \c BoxGenerator
instance) and query boxes of the same size (second \c BoxGenerator instance).

\c libfbi now generates an adjacency list \c adjList that describes a graph in
which all centroids (the vertices of the graph) whose uncertainty boxes overlap
are connected with an edge. Determining the XICs now is a simple
matter of determining the connected components of the graph:

\code
    typedef SetA<Centroid, 1, 2>::IntType LabelType;
    std::vector<LabelType> labels;
    findConnectedComponents(adjList, labels); 
\endcode

Finally, we can write the results to a text file and exit:
\code
    std::ofstream ofs(options.outputfileName_.c_str());
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    for (size_t i = 0;i < centroids.size();++i) {
    Centroid& c = centroids[i];
    ofs << c.rt_ << "\t" << c.mz_ << "\t" << c.sn_ << "\t" 
      << c.abundance_ << "\t" << labels[i] << "\n";
    }
    return 0;
}
\endcode
*/
