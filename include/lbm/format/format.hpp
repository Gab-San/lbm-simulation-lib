#pragma once

/// Minimal brace-substitution formatting over std::ostream. No dependencies.
///
/// Header-only on purpose: the ostream log backend and lbm::csv both use it,
/// and neither can gain a link dependency without editing the CMake module
/// that builds it. Everything here is inline, so including it is enough.
///
/// Supported:
///   {}            positional, in order
///   {{  }}        literal braces
///   {:.3f}        precision + float style (f F e E g G)
///   {:x} {:#X}    integer base (x X o d), '#' for showbase
///   {:+}          showpos
///
/// Parsed and deliberately ignored (never an error): fill, alignment, width,
/// zero-padding, and any unrecognised type character.
///
/// Malformed input degrades to literal output rather than throwing: an
/// unterminated '{' emits the remainder of the format string verbatim, and a
/// placeholder with no corresponding argument expands to nothing.
///
/// This header must not include anything from lbm-sim/.

#include <cstddef>
#include <ios>
#include <memory>
#include <ostream>
#include <string_view>

namespace lbm::format {

namespace detail {

/// A type-erased argument: a pointer to the caller's object plus the function
/// that knows how to print it. No allocation, no std::function.
struct Arg {
  void const *value;
  void (*print)(std::ostream &, void const *, std::string_view);
};

/// Saves and restores stream formatting state around a single argument, so a
/// spec on one placeholder cannot leak into the next one or into the caller's
/// stream.
class StateGuard {
public:
  explicit StateGuard(std::ostream &os)
      : os_(os), flags_(os.flags()), precision_(os.precision()) {}

  ~StateGuard() {
    os_.flags(flags_);
    os_.precision(precision_);
  }

  StateGuard(StateGuard const &) = delete;
  StateGuard &operator=(StateGuard const &) = delete;

private:
  std::ostream &os_;
  std::ios_base::fmtflags flags_;
  std::streamsize precision_;
};

constexpr bool is_digit(char c) noexcept { return c >= '0' && c <= '9'; }

constexpr bool is_align(char c) noexcept {
  return c == '<' || c == '>' || c == '^';
}

/// Translates the text between ':' and '}' into stream manipulator calls.
inline void apply_spec(std::ostream &os, std::string_view spec) {
  std::size_t i = 0;

  // [[fill]align] - consumed, not honoured
  if (i + 1 < spec.size() && is_align(spec[i + 1])) {
    i += 2;
  } else if (i < spec.size() && is_align(spec[i])) {
    i += 1;
  }

  // [sign]
  if (i < spec.size()) {
    if (spec[i] == '+') {
      os << std::showpos;
      ++i;
    } else if (spec[i] == '-' || spec[i] == ' ') {
      ++i;
    }
  }

  // [#]
  if (i < spec.size() && spec[i] == '#') {
    os << std::showbase;
    ++i;
  }

  // [0][width] - consumed, not honoured
  while (i < spec.size() && is_digit(spec[i])) {
    ++i;
  }

  // [.precision]
  if (i < spec.size() && spec[i] == '.') {
    ++i;
    std::streamsize precision = 0;
    bool seen = false;
    while (i < spec.size() && is_digit(spec[i])) {
      precision = precision * 10 + (spec[i] - '0');
      seen = true;
      ++i;
    }
    if (seen) {
      os.precision(precision);
    }
  }

  // [type]
  if (i < spec.size()) {
    switch (spec[i]) {
    case 'f':
    case 'F':
      os << std::fixed;
      break;
    case 'e':
      os << std::scientific << std::nouppercase;
      break;
    case 'E':
      os << std::scientific << std::uppercase;
      break;
    case 'g':
      os.unsetf(std::ios_base::floatfield);
      os << std::nouppercase;
      break;
    case 'G':
      os.unsetf(std::ios_base::floatfield);
      os << std::uppercase;
      break;
    case 'x':
      os << std::hex << std::nouppercase;
      break;
    case 'X':
      os << std::hex << std::uppercase;
      break;
    case 'o':
      os << std::oct;
      break;
    case 'd':
      os << std::dec;
      break;
    default:
      break;
    }
  }
}

template <class T>
void print_arg(std::ostream &os, void const *value, std::string_view spec) {
  auto const &v = *static_cast<T const *>(value);
  if (spec.empty()) {
    os << v;
    return;
  }
  StateGuard const guard(os);
  apply_spec(os, spec);
  os << v;
}

template <class T> Arg make_arg(T const &v) noexcept {
  return Arg{static_cast<void const *>(std::addressof(v)), &print_arg<T>};
}

inline void vformat_to(std::ostream &os, std::string_view fmt, Arg const *args,
                       std::size_t count) {
  constexpr auto npos = std::string_view::npos;

  std::size_t next = 0;
  std::size_t i = 0;

  while (i < fmt.size()) {
    std::size_t const brace = fmt.find_first_of("{}", i);
    if (brace == npos) {
      os << fmt.substr(i);
      return;
    }
    if (brace > i) {
      os << fmt.substr(i, brace - i);
    }

    char const c = fmt[brace];

    // "{{" or "}}"
    if (brace + 1 < fmt.size() && fmt[brace + 1] == c) {
      os << c;
      i = brace + 2;
      continue;
    }

    // Stray '}' - emit literally.
    if (c == '}') {
      os << c;
      i = brace + 1;
      continue;
    }

    std::size_t const close = fmt.find('}', brace + 1);
    if (close == npos) {
      os << fmt.substr(brace);
      return;
    }

    std::string_view const body = fmt.substr(brace + 1, close - brace - 1);
    std::size_t const colon = body.find(':');
    std::string_view const spec =
        (colon == npos) ? std::string_view{} : body.substr(colon + 1);

    if (next < count) {
      Arg const &arg = args[next++];
      arg.print(os, arg.value, spec);
    }

    i = close + 1;
  }
}

} // namespace detail

/// Writes `fmt` to `os`, substituting `args` for successive placeholders.
/// Requires only that each argument has an `operator<<` for std::ostream.
template <class... Ts>
void format_to(std::ostream &os, std::string_view fmt, Ts const &...args) {
  if constexpr (sizeof...(Ts) == 0) {
    detail::vformat_to(os, fmt, nullptr, 0);
  } else {
    detail::Arg const argv[] = {detail::make_arg(args)...};
    detail::vformat_to(os, fmt, argv, sizeof...(Ts));
  }
}

} // namespace lbm::format
