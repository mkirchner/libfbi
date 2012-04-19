#include <tuple>
template <typename ... T>
struct test : public std::tuple<T...>
{};
int main() {
test<int, int, int> a;

}
