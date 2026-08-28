#ifndef MAGNETIC_COORDINATE_DETAIL_FIND_ROOT_HPP_
#define MAGNETIC_COORDINATE_DETAIL_FIND_ROOT_HPP_

// TOMS748 implementation adapted from libmeq util.h by Zhong Qi, distributed
// under the MIT license in the repository root. The interpolation sequence is
// intentionally retained so magnetic-coordinate and libmeq choose the same
// roots for the same bracket and tolerance.

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace magnetic_coordinate::detail {

template <typename T>
constexpr T root_abs(T value) {
    return value < T{0} ? -value : value;
}

template <typename T>
constexpr int root_sign(const T& value) {
    return (T{0} < value) - (value < T{0});
}

template <typename Function, typename X, typename F>
void bracket(Function& function,
             X& left,
             X& right,
             X candidate,
             F& f_left,
             F& f_right,
             X& previous,
             F& f_previous) {
    const F tolerance = std::numeric_limits<F>::epsilon() * F{2};
    if ((right - left) < X{2} * tolerance * left) {
        candidate = left + (right - left) / X{2};
    } else if (candidate - left <= root_abs(left) * tolerance) {
        candidate = left + root_abs(left) * tolerance;
    } else if (right - candidate <= root_abs(right) * tolerance) {
        candidate = right - root_abs(right) * tolerance;
    }

    const F f_candidate = function(candidate);
    if (std::fpclassify(f_candidate) == FP_ZERO) {
        left = candidate;
        f_left = F{0};
        previous = X{0};
        f_previous = F{0};
        return;
    }
    if (root_sign(f_left) * root_sign(f_candidate) < 0) {
        previous = right;
        f_previous = f_right;
        right = candidate;
        f_right = f_candidate;
    } else {
        previous = left;
        f_previous = f_left;
        left = candidate;
        f_left = f_candidate;
    }
}

template <typename T>
T safe_divide(T numerator, T denominator, T fallback) {
    return root_abs(denominator) < T{1} &&
                   root_abs(denominator * std::numeric_limits<T>::max()) <=
                       root_abs(numerator)
               ? fallback
               : numerator / denominator;
}

template <typename X, typename F>
X secant_interpolate(const X& left,
                     const X& right,
                     const F& f_left,
                     const F& f_right) {
    const F tolerance = std::numeric_limits<F>::epsilon() * F{5};
    const X candidate = left - (f_left / (f_right - f_left)) * (right - left);
    return candidate - left <= root_abs(left) * tolerance ||
                   right - candidate <= root_abs(right) * tolerance
               ? (left + right) / X{2}
               : candidate;
}

template <typename T>
bool outside(const T& candidate, const T& left, const T& right) {
    return candidate <= left || candidate >= right;
}

template <typename X, typename F>
X quadratic_interpolation(const X& left,
                          const X& right,
                          const X& previous,
                          const F& f_left,
                          const F& f_right,
                          const F& f_previous,
                          unsigned count) {
    F linear = safe_divide(F(f_right - f_left), F(right - left),
                           std::numeric_limits<F>::max());
    F quadratic = safe_divide(F(f_previous - f_right), F(previous - right),
                              std::numeric_limits<F>::max());
    quadratic = safe_divide(F(quadratic - linear), F(previous - left), F{0});
    if (std::fpclassify(quadratic) == FP_ZERO) {
        return secant_interpolate(left, right, f_left, f_right);
    }

    X candidate = root_sign(quadratic) * root_sign(f_left) > 0 ? left : right;
    for (unsigned iteration = 1; iteration <= count; ++iteration) {
        candidate -= safe_divide(
            F(f_left +
              (linear + quadratic * (candidate - right)) * (candidate - left)),
            F(linear + quadratic * (X{2} * candidate - left - right)),
            F{1} + candidate - left);
    }
    return outside(candidate, left, right)
               ? secant_interpolate(left, right, f_left, f_right)
               : candidate;
}

template <typename X, typename F>
X cubic_interpolation(const X& left,
                      const X& right,
                      const X& previous,
                      const X& older,
                      const F& f_left,
                      const F& f_right,
                      const F& f_previous,
                      const F& f_older) {
    const X q11 = (previous - older) * f_previous / (f_older - f_previous);
    const X q21 = (right - previous) * f_right / (f_previous - f_right);
    const X q31 = (left - right) * f_left / (f_right - f_left);
    const X d21 = (right - previous) * f_previous / (f_previous - f_right);
    const X d31 = (left - right) * f_right / (f_right - f_left);
    const X q22 = (d21 - q11) * f_right / (f_older - f_right);
    const X q32 = (d31 - q21) * f_left / (f_previous - f_left);
    const X d32 = (d31 - q21) * f_previous / (f_previous - f_left);
    const X q33 = (d32 - q22) * f_left / (f_older - f_left);
    const X candidate = q31 + q32 + q33 + left;
    return outside(candidate, left, right)
               ? quadratic_interpolation(left, right, previous, f_left, f_right,
                                         f_previous, 3)
               : candidate;
}

template <typename Function, typename Tolerance, typename X>
X find_root(const Function& function,
            const X& left_bound,
            const X& right_bound,
            const Tolerance& tolerance) {
    using F = std::invoke_result_t<Function, X>;
    constexpr unsigned MAX_ITERATIONS = 50;
    unsigned remaining = MAX_ITERATIONS;
    constexpr F MINIMUM_CONTRACTION{0.5L};

    if (!(left_bound < right_bound)) {
        throw std::domain_error("root interval does not exist");
    }
    X left = left_bound;
    X right = right_bound;
    F f_left = function(left);
    F f_right = function(right);
    if (!std::isfinite(f_left) || !std::isfinite(f_right) ||
        root_sign(f_left) * root_sign(f_right) > 0) {
        throw std::domain_error("root interval does not bracket a finite root");
    }
    if (tolerance(left, right) || std::fpclassify(f_left) == FP_ZERO ||
        std::fpclassify(f_right) == FP_ZERO) {
        return std::fpclassify(f_right) == FP_ZERO ? right : left;
    }

    X previous{0};
    X older{1.0e5L};
    F f_previous{1.0e5L};
    F f_older{1.0e5L};
    X candidate = secant_interpolate(left, right, f_left, f_right);
    auto mutable_function = function;
    bracket(mutable_function, left, right, candidate, f_left, f_right, previous,
            f_previous);
    --remaining;
    if (remaining != 0 && std::fpclassify(f_left) != FP_ZERO &&
        !tolerance(left, right)) {
        candidate = quadratic_interpolation(left, right, previous, f_left,
                                            f_right, f_previous, 2);
        older = previous;
        f_older = f_previous;
        bracket(mutable_function, left, right, candidate, f_left, f_right,
                previous, f_previous);
        --remaining;
    }

    const F minimum_difference = std::numeric_limits<F>::min() * F{32};
    while (remaining != 0 && std::fpclassify(f_left) != FP_ZERO &&
           !tolerance(left, right)) {
        const X old_left = left;
        const X old_right = right;
        auto values_too_close = [&]() {
            return root_abs(f_left - f_right) < minimum_difference ||
                   root_abs(f_left - f_previous) < minimum_difference ||
                   root_abs(f_left - f_older) < minimum_difference ||
                   root_abs(f_right - f_previous) < minimum_difference ||
                   root_abs(f_right - f_older) < minimum_difference ||
                   root_abs(f_previous - f_older) < minimum_difference;
        };

        candidate =
            values_too_close()
                ? quadratic_interpolation(left, right, previous, f_left,
                                          f_right, f_previous, 2)
                : cubic_interpolation(left, right, previous, older, f_left,
                                      f_right, f_previous, f_older);
        older = previous;
        f_older = f_previous;
        bracket(mutable_function, left, right, candidate, f_left, f_right,
                previous, f_previous);
        if (--remaining == 0 || std::fpclassify(f_left) == FP_ZERO ||
            tolerance(left, right)) {
            break;
        }

        candidate =
            values_too_close()
                ? quadratic_interpolation(left, right, previous, f_left,
                                          f_right, f_previous, 3)
                : cubic_interpolation(left, right, previous, older, f_left,
                                      f_right, f_previous, f_older);
        bracket(mutable_function, left, right, candidate, f_left, f_right,
                previous, f_previous);
        if (--remaining == 0 || std::fpclassify(f_left) == FP_ZERO ||
            tolerance(left, right)) {
            break;
        }

        X best = right;
        F f_best = f_right;
        if (root_abs(f_left) < root_abs(f_right)) {
            best = left;
            f_best = f_left;
        }
        candidate =
            best - X{2} * (f_best / (f_right - f_left)) * (right - left);
        if (root_abs(candidate - best) > (right - left) / X{2}) {
            candidate = left + (right - left) / X{2};
        }
        older = previous;
        f_older = f_previous;
        bracket(mutable_function, left, right, candidate, f_left, f_right,
                previous, f_previous);
        if (--remaining == 0 || std::fpclassify(f_left) == FP_ZERO ||
            tolerance(left, right)) {
            break;
        }
        if ((right - left) < MINIMUM_CONTRACTION * (old_right - old_left)) {
            continue;
        }
        older = previous;
        f_older = f_previous;
        bracket(mutable_function, left, right, left + (right - left) / X{2},
                f_left, f_right, previous, f_previous);
        --remaining;
    }
    return std::fpclassify(f_right) == FP_ZERO ? right : left;
}

}  // namespace magnetic_coordinate::detail

#endif  // MAGNETIC_COORDINATE_DETAIL_FIND_ROOT_HPP_
