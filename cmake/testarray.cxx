#include <array>

int main() {
  typedef std::array<int, 3> TupleType;
  TupleType a;
  std::get<0>(a);
  std::tuple_size<TupleType>::value;
}
