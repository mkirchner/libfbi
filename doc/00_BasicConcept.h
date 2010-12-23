/*!
\page basic_concept Basic Concepts
The main functionality of \c libfbi is the fast determination of pairs of
intersection boxes from two sets. With respect to the interface, this
functionality is provided by a single static function:
\code
fbi::setA<BoxType1, TIndices1>::setB<BoxType2, TIndices2>::intersect(a, aBoxGen,  b, bBoxGen)
\endcode
where \c a is the container (e.g. \c std::vector) of elements that make up the
first set and \c b is the container that makes up the second set. The function
returns a \c std::list of \c std::sets that holds the adjecency information for
each element in both sets (the adjacency list index starts at 0, indexes all
elements in \c a, then continues indexing all elements in \c b).

\section sec_intersect Using intersect(...)
It is important to realize that \c libfbi never touches the data contained in
\c a or \c b. Instead, it constructs a box for each element in \c a and \c b on
the fly and uses these boxes throughout the calculation. But how does \c libfbi
acquire this information? This is where the user comes in: it is the task of
the user to provide \c libfbi with functions that deliver the interval
information in each dimension. Assume we have a type \c X. \c X has many data
members, but for the context of the range search only a few of them are
relevant:
\code
struct X {
    double a, b, c, d, e, f, g, h, i, j;
    double sigma_b, sigma_h;
}
\endcode
Assuming that the box intersection should be conducted for \c b and \c h, the
user would provide the helper functions that extract the relevant information:
\code
class BoxGeneratorX
{
  public:
    BoxGeneratorX(...) : ... {}

    template <size_t N>
      typename std::tuple_element<N,
      typename fbi::Traits<X>::key_type>::type
    get(const X&) const;
};

template <>
std::pair<double, double>
BoxGeneratorX::get<0>(const X& x) const
{
    return std::make_pair(x.b-3*x.sigma_b, x.b+3*x.sigma_b);
}

template <>
std::pair<double, double>
BoxGeneratorX::get<1>(const X& x) const
{
    return std::make_pair(x.g-3*x.sigma_g, x.g+3*x.sigma_g);
}
\endcode
The \c get<0> \a accessor function now returns an interval for \c x.b and \c
get<1> returns an interval for \c x.h. Assuming that we would like to intersect
elements of type \c X with elements of type \c Y, we also need a definition for
\c Y:
\code
struct Y {
    double m, n, o, p, q;
};
\endcode
The accessor functions for \c BoxGeneratorY can now be implemented similar to
those for \c X, e.g:
\code
class BoxGeneratorY
{
  public:
    BoxGeneratorY(...) : ... {}

    template <size_t N>
      typename std::tuple_element<N,
      typename fbi::Traits<Y>::key_type>::type
    get(const Y&) const;
};

template <>
std::pair<double, double>
BoxGeneratorY::get<0>(const Y& y) const
{
    return std::make_pair(y.m-7.54, y.m+1.3);
}

template <>
std::pair<double, double>
BoxGeneratorY::get<1>(const Y& y) const
{
    // this interval is completely to the left of y.p; this is possible.
    return std::make_pair(y.p-3.14, y.p-1.14);
}
\endcode
Once the respective declarations and definitions are in place,
calling the fast box intersection procedure is simply accomplished by
\code
auto adjList = fbi::setA<X, 1, 2>::setB<Y, 1, 2>::intersect(x, BoxGeneratorX,
y, BoxGeneratorY);
\endcode
where \c x and \c y are containers of instances of \c X and \c Y, respectively.

In the common use case where one is interested in self-intersection of a set of
boxes, \c libfbi provides a useful shortcut. Assume one wants to determine the
set of boxes in \c a that overlap with each other, given the box definition for
the \c X type from above. In that case, the expression simplifies to
\code
auto adjList = fbi::setA<X, 1, 2>::intersect(x, BoxGeneratorX);
\endcode

\section sec_adjlist Using the adjacency list
The call to \c intersect(...) returns an adjacency list, which in general is a
\c std::vector<std::set<size_t> > object. The adjacency list has size \c
x.size() + \c y.size(), is indexed starting from 0, runs through the
indexes of all elements of \c x, and continues to index all elements of \c y.
Hence the query results for \c x[k] are given by
\code
...
std::set<size_t> queryResultIndexes = adjList[k];
...
\endcode
However, note that the indexes for all elements in \c y are shifted by \c
x.size(). Hence, indexing the elements of \c y requires some index
calculation. For \c y[k], the query result is available with:
\code
...
std::set<size_t> queryResultIndexes = adjList[k+x.size()];
...
\endcode
The following illustrates indexing elements in \c y given its index in the 
adjacency list:
\code
...
size_t shift = x.size();
typedef std::vector<std::set<size_t> >::const_iterator VSCI;
for (VSCI i = adjList.begin(); i != adjList.end(); ++i) {
    std::cout << x[std::distance(adjList.begin(), i)] << '\n';
    typedef std::set<size_t>::const_iterator SCI; 
    for (SCI j = i->begin(); j != i->end(); ++j) {
        std::cout << "  +- " << y[*j - shift] << '\n';
    }
}
\endcode
This code will print all elements in \c x and \c y, each with all its neighbors
indented beneath it. Note that (because box intersection calculates undirected
graphs) this prints each pair twice: for every
combination (m,n) present in \c adjList, there is also the combination (n,m).
*/
