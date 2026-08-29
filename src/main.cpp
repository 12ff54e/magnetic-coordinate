#include "magnetic_coordinate/boozer_output.hpp"
#include "magnetic_coordinate/transform.hpp"

#include <charconv>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

struct CommandLine {
    std::filesystem::path input;
    std::filesystem::path output = "boozer-output.bin";
    magnetic_coordinate::TransformSettings settings;
};

int parse_integer(std::string_view text, const char* option) {
    int value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        throw std::invalid_argument(std::string(option) +
                                    " requires an integer");
    }
    return value;
}

double parse_double(std::string_view text, const char* option) {
    double value = 0.0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        throw std::invalid_argument(std::string(option) + " requires a number");
    }
    return value;
}

CommandLine parse_command_line(int argc, char** argv) {
    CommandLine command;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string_view token(argv[argument]);
        const auto value_after = [&](const char* option) -> std::string_view {
            if (++argument >= argc) {
                throw std::invalid_argument(std::string(option) +
                                            " requires a value");
            }
            return argv[argument];
        };
        if (token == "-o" || token == "--output") {
            command.output = value_after("--output");
        } else if (token == "--ntheta") {
            command.settings.output_ntheta =
                parse_integer(value_after("--ntheta"), "--ntheta");
        } else if (token == "--nzeta") {
            command.settings.output_nzeta =
                parse_integer(value_after("--nzeta"), "--nzeta");
        } else if (token == "--mmax") {
            command.settings.mmax =
                parse_integer(value_after("--mmax"), "--mmax");
        } else if (token == "--nmax") {
            command.settings.nmax =
                parse_integer(value_after("--nmax"), "--nmax");
        } else if (token == "--radial-order") {
            const int order =
                parse_integer(value_after("--radial-order"), "--radial-order");
            if (order == 2) {
                command.settings.radial_order =
                    magnetic_coordinate::RadialInterpolationOrder::TWO_POINT;
            } else if (order == 4) {
                command.settings.radial_order =
                    magnetic_coordinate::RadialInterpolationOrder::FOUR_POINT;
            } else {
                throw std::invalid_argument("--radial-order must be 2 or 4");
            }
        } else if (token == "--resonance-tolerance") {
            command.settings.resonance_tolerance = parse_double(
                value_after("--resonance-tolerance"), "--resonance-tolerance");
        } else if (token == "-h" || token == "--help") {
            std::cout
                << "usage: cumes-boozer INPUT [--output FILE] [--ntheta N] "
                   "[--nzeta N] "
                   "[--mmax M] [--nmax N] [--radial-order 2|4] "
                   "[--resonance-tolerance X]\n";
            std::exit(0);
        } else if (!token.empty() && token.front() == '-') {
            throw std::invalid_argument("unknown option: " +
                                        std::string(token));
        } else if (command.input.empty()) {
            command.input = token;
        } else {
            throw std::invalid_argument("only one input file is accepted");
        }
    }
    if (command.input.empty()) {
        throw std::invalid_argument("a cuMES input file is required");
    }
    return command;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto command = parse_command_line(argc, argv);
        const auto output_spec =
            magnetic_coordinate::resolve_boozer_output_spec(command.output);
        if (!magnetic_coordinate::boozer_output_format_available(
                output_spec.format)) {
            throw std::runtime_error(
                "output format '" +
                std::string(magnetic_coordinate::boozer_output_suffix(
                    output_spec.format)) +
                "' is not available in this build");
        }
        const auto result = magnetic_coordinate::transform_cumes_file(
            command.input, command.settings);
        magnetic_coordinate::write_boozer_output(output_spec, result,
                                                 command.input);
        const auto real_modes =
            static_cast<std::size_t>(result.spectrum.mmax + 1) *
            static_cast<std::size_t>(result.spectrum.nmax + 1);
        std::cout << "wrote " << command.output << " with " << result.s.size()
                  << " non-axis surfaces, " << result.grid.ntheta << "x"
                  << result.grid.nzeta << " mixed-grid points per surface, and "
                  << real_modes << " real modes per surface\n";
    } catch (const std::exception& error) {
        std::cerr << "cumes-boozer: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
