#include <tuple>
template <typename ...T>
using test = std::tupl<T...>;
int main() {
test<int, int, int> a;
}
