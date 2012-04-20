#include <array>
struct kdtree {};
template <typename Tree>
struct kdtreeTraits {
typedef Tree tt;
};
template <typename Traits>
struct Dummy {
friend class Traits::tt;

};


int main() {
  typedef std::array<int, 3> TupleType;
  TupleType{{2, 3, 4}};
  TupleType a;
  std::get<0>(a);
  std::tuple_size<TupleType>::value;
  Dummy<kdtreeTraits<kdtree> > dummy;
}
