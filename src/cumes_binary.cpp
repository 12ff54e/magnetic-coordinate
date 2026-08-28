#include "magnetic_coordinate/cumes_binary.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace magnetic_coordinate {
namespace {

constexpr std::array<char, 8> MAGIC{'C', 'U', 'M', 'E', 'S', '0', '0', '1'};
constexpr std::int32_t REQUIRED_VERSION = 8;
constexpr std::int32_t MAX_VECTOR_ELEMENTS = 1 << 20;
constexpr std::size_t MAX_REAL_FIELD_ELEMENTS = 1U << 24;
constexpr std::int32_t MAX_STRING_BYTES = 1 << 24;

std::size_t checked_product(std::size_t left,
                            std::size_t right,
                            const char* description) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::runtime_error(std::string("cuMES binary: ") + description +
                                 " overflows size_t");
    }
    return left * right;
}

class BinaryReader {
   public:
    explicit BinaryReader(const std::filesystem::path& path)
        : input_(path, std::ios::binary), path_(path.string()) {
        if (!input_) {
            throw std::runtime_error("cannot open cuMES binary " + path_);
        }
    }

    std::array<char, 8> read_magic() {
        std::array<char, 8> value{};
        read_exact(value.data(), value.size(), "magic");
        return value;
    }

    std::uint8_t read_u8(const char* description) {
        std::uint8_t value = 0;
        read_exact(&value, sizeof(value), description);
        return value;
    }

    std::int32_t read_i32(const char* description) {
        const std::uint32_t bits = read_little_uint<std::uint32_t>(description);
        return std::bit_cast<std::int32_t>(bits);
    }

    double read_f64(const char* description) {
        const std::uint64_t bits = read_little_uint<std::uint64_t>(description);
        return std::bit_cast<double>(bits);
    }

    std::string read_string(const char* description) {
        const std::int32_t count = read_i32(description);
        if (count < 0 || count > MAX_STRING_BYTES) {
            throw std::runtime_error(std::string("cuMES binary: invalid ") +
                                     description + " length");
        }
        std::string value(static_cast<std::size_t>(count), '\0');
        read_exact(value.data(), value.size(), description);
        return value;
    }

    std::vector<double> read_f64_vector(const char* description) {
        const std::int32_t count = read_bounded_count(description);
        std::vector<double> values(static_cast<std::size_t>(count));
        for (double& value : values) value = read_f64(description);
        return values;
    }

    void skip_i32_vector(const char* description) {
        const std::int32_t count = read_bounded_count(description);
        skip_bytes(checked_product(static_cast<std::size_t>(count),
                                   sizeof(std::int32_t), description),
                   description);
    }

    void skip_f64_vector(const char* description) {
        const std::int32_t count = read_bounded_count(description);
        skip_bytes(checked_product(static_cast<std::size_t>(count),
                                   sizeof(double), description),
                   description);
    }

    std::vector<double> read_f64_array(std::size_t count,
                                       const char* description) {
        std::vector<double> values(count);
        for (double& value : values) value = read_f64(description);
        return values;
    }

    void skip_bytes(std::size_t count, const char* description) {
        if (count > static_cast<std::size_t>(
                        std::numeric_limits<std::streamoff>::max())) {
            throw std::runtime_error(std::string("cuMES binary: ") +
                                     description + " is too large");
        }
        input_.seekg(static_cast<std::streamoff>(count), std::ios::cur);
        if (!input_) truncated(description);
    }

   private:
    template <typename UInt>
    UInt read_little_uint(const char* description) {
        static_assert(std::is_unsigned_v<UInt>);
        std::array<std::uint8_t, sizeof(UInt)> bytes{};
        read_exact(bytes.data(), bytes.size(), description);
        UInt value = 0;
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            value |= static_cast<UInt>(bytes[i]) << (8U * i);
        }
        return value;
    }

    std::int32_t read_bounded_count(const char* description) {
        const std::int32_t count = read_i32(description);
        if (count < 0 || count > MAX_VECTOR_ELEMENTS) {
            throw std::runtime_error(std::string("cuMES binary: invalid ") +
                                     description + " count");
        }
        return count;
    }

    void read_exact(void* destination,
                    std::size_t count,
                    const char* description) {
        if (count == 0) return;
        input_.read(static_cast<char*>(destination),
                    static_cast<std::streamsize>(count));
        if (!input_) truncated(description);
    }

    [[noreturn]] void truncated(const char* description) const {
        throw std::runtime_error("cuMES binary: truncated " +
                                 std::string(description) + " in " + path_);
    }

    std::ifstream input_;
    std::string path_;
};

void skip_output_provenance(BinaryReader& reader) {
    static_cast<void>(reader.read_i32("precision tag"));
    static_cast<void>(reader.read_i32("run status"));
    static_cast<void>(reader.read_i32("total iteration count"));
    const std::int32_t stage_count = reader.read_i32("stage count");
    if (stage_count < 0 || stage_count > MAX_VECTOR_ELEMENTS) {
        throw std::runtime_error("cuMES binary: invalid stage count");
    }

    static_cast<void>(reader.read_string("build revision"));
    static_cast<void>(reader.read_u8("dirty flag"));
    static_cast<void>(reader.read_string("build type"));
    static_cast<void>(reader.read_string("precision policy"));
    static_cast<void>(reader.read_string("compile flags"));
    static_cast<void>(reader.read_string("input path"));
    static_cast<void>(reader.read_string("input hash"));
    static_cast<void>(reader.read_string("GPU name"));
    static_cast<void>(reader.read_string("driver version"));
    static_cast<void>(reader.read_string("runtime version"));
    static_cast<void>(reader.read_string("toolkit version"));

    for (std::int32_t stage = 0; stage < stage_count; ++stage) {
        static_cast<void>(reader.read_i32("stage ns"));
        static_cast<void>(reader.read_i32("stage iterations"));
        static_cast<void>(reader.read_u8("stage converged flag"));
        static_cast<void>(reader.read_f64("stage fsqr"));
        static_cast<void>(reader.read_f64("stage fsqz"));
        static_cast<void>(reader.read_f64("stage fsql"));
        const std::int32_t restarts = reader.read_i32("restart count");
        if (restarts < 0 || restarts > MAX_VECTOR_ELEMENTS) {
            throw std::runtime_error("cuMES binary: invalid restart count");
        }
        reader.skip_bytes(checked_product(static_cast<std::size_t>(restarts),
                                          sizeof(std::int32_t), "restart list"),
                          "restart list");
    }
}

void read_input_params(BinaryReader& reader, CumesEquilibrium& equilibrium) {
    equilibrium.mpol = reader.read_i32("mpol");
    equilibrium.ntor = reader.read_i32("ntor");
    equilibrium.nfp = reader.read_i32("nfp");
    const std::int32_t input_ntheta = reader.read_i32("input ntheta");
    const std::int32_t input_nzeta = reader.read_i32("input nzeta");
    equilibrium.ncurr = reader.read_i32("ncurr");

    static_cast<void>(reader.read_f64("delt"));
    equilibrium.phiedge = reader.read_f64("phiedge");
    for (int scalar = 0; scalar < 6; ++scalar) {
        static_cast<void>(reader.read_f64("input scalar"));
    }

    static_cast<void>(reader.read_string("input schema"));
    static_cast<void>(reader.read_string("pressure profile type"));
    static_cast<void>(reader.read_string("iota profile type"));
    static_cast<void>(reader.read_string("current profile type"));

    reader.skip_f64_vector("pressure coefficients");
    reader.skip_f64_vector("current coefficients");
    reader.skip_f64_vector("iota coefficients");
    equilibrium.aphi = reader.read_f64_vector("toroidal flux coefficients");
    reader.skip_f64_vector("R-axis coefficients");
    reader.skip_f64_vector("Z-axis coefficients");

    const std::int32_t input_stages = reader.read_i32("input stage count");
    if (input_stages < 0 || input_stages > MAX_VECTOR_ELEMENTS) {
        throw std::runtime_error("cuMES binary: invalid input stage count");
    }
    for (std::int32_t stage = 0; stage < input_stages; ++stage) {
        static_cast<void>(reader.read_i32("input stage ns"));
        static_cast<void>(reader.read_i32("input stage iterations"));
        static_cast<void>(reader.read_f64("input stage tolerance"));
    }

    reader.skip_i32_vector("rbc m");
    reader.skip_i32_vector("rbc n");
    reader.skip_f64_vector("rbc values");
    reader.skip_i32_vector("zbs m");
    reader.skip_i32_vector("zbs n");
    reader.skip_f64_vector("zbs values");
    reader.skip_f64_vector("rbcc");
    reader.skip_f64_vector("rbss");
    reader.skip_f64_vector("zbsc");
    reader.skip_f64_vector("zbcs");

    static_cast<void>(reader.read_i32("free-boundary flag"));
    static_cast<void>(reader.read_i32("vacuum skip"));
    static_cast<void>(reader.read_string("mgrid path"));
    reader.skip_f64_vector("external currents");
    static_cast<void>(reader.read_string("coils path"));
    static_cast<void>(reader.read_string("Makegrid parameters path"));
    const std::int32_t embedded = reader.read_i32("embedded Makegrid flag");
    if (embedded != 0) {
        static_cast<void>(reader.read_i32("Makegrid normalize flag"));
        static_cast<void>(reader.read_i32("Makegrid symmetry flag"));
        static_cast<void>(reader.read_i32("Makegrid field periods"));
        static_cast<void>(reader.read_f64("Makegrid R minimum"));
        static_cast<void>(reader.read_f64("Makegrid R maximum"));
        static_cast<void>(reader.read_i32("Makegrid R points"));
        static_cast<void>(reader.read_f64("Makegrid Z minimum"));
        static_cast<void>(reader.read_f64("Makegrid Z maximum"));
        static_cast<void>(reader.read_i32("Makegrid Z points"));
        static_cast<void>(reader.read_i32("Makegrid phi points"));
    }

    if (equilibrium.mpol < 1 || equilibrium.ntor < 0 || equilibrium.nfp < 1 ||
        input_ntheta != equilibrium.ntheta ||
        input_nzeta != equilibrium.nzeta || equilibrium.aphi.empty()) {
        throw std::runtime_error(
            "cuMES binary: inconsistent or incomplete input parameters");
    }
    const std::size_t expected_modes = checked_product(
        static_cast<std::size_t>(equilibrium.mpol),
        static_cast<std::size_t>(equilibrium.ntor + 1), "mode count");
    if (expected_modes != static_cast<std::size_t>(equilibrium.mnmax)) {
        throw std::runtime_error(
            "cuMES binary: mnmax disagrees with mpol and ntor");
    }
}

}  // namespace

NativeFieldView CumesEquilibrium::native_field_view() const {
    return NativeFieldView{
        ns,
        ntheta,
        nzeta,
        half_fields[SQRTG],
        half_fields[BSUPS],
        half_fields[BSUPU],
        half_fields[BSUPV],
        half_fields[BSUBS],
        half_fields[BSUBU],
        half_fields[BSUBV],
    };
}

CumesEquilibriumView CumesEquilibrium::view() const {
    CumesEquilibriumView result;
    result.format_version = format_version;
    result.ns = ns;
    result.mnmax = mnmax;
    result.mpol = mpol;
    result.ntor = ntor;
    result.nfp = nfp;
    result.ntheta = ntheta;
    result.nzeta = nzeta;
    result.ncurr = ncurr;
    result.phiedge = phiedge;
    result.aphi = aphi;
    for (std::size_t family = 0; family < families.size(); ++family) {
        result.families[family] = families[family];
    }
    for (std::size_t field = 0; field < half_fields.size(); ++field) {
        result.half_fields[field] = half_fields[field];
    }
    return result;
}

NativeFieldView CumesEquilibriumView::native_field_view() const {
    return NativeFieldView{
        ns,
        ntheta,
        nzeta,
        half_fields[CumesEquilibrium::SQRTG],
        half_fields[CumesEquilibrium::BSUPS],
        half_fields[CumesEquilibrium::BSUPU],
        half_fields[CumesEquilibrium::BSUPV],
        half_fields[CumesEquilibrium::BSUBS],
        half_fields[CumesEquilibrium::BSUBU],
        half_fields[CumesEquilibrium::BSUBV],
    };
}

CumesEquilibrium read_cumes_binary(const std::filesystem::path& path) {
    BinaryReader reader(path);
    if (reader.read_magic() != MAGIC) {
        throw std::runtime_error("cuMES binary: bad magic");
    }

    CumesEquilibrium equilibrium;
    equilibrium.format_version = reader.read_i32("format version");
    if (equilibrium.format_version != REQUIRED_VERSION) {
        throw std::runtime_error(
            "cuMES binary: magnetic-coordinate requires format version 8");
    }
    equilibrium.ns = reader.read_i32("ns");
    equilibrium.mnmax = reader.read_i32("mnmax");
    if (equilibrium.ns < 3 || equilibrium.mnmax < 1) {
        throw std::runtime_error("cuMES binary: invalid spectral dimensions");
    }

    const std::size_t family_size = checked_product(
        static_cast<std::size_t>(equilibrium.ns),
        static_cast<std::size_t>(equilibrium.mnmax), "spectral dimensions");
    for (auto& family : equilibrium.families) {
        family = reader.read_f64_array(family_size, "spectral family");
    }

    equilibrium.ntheta = reader.read_i32("ntheta");
    equilibrium.nzeta = reader.read_i32("nzeta");
    if (equilibrium.ntheta < 1 || equilibrium.nzeta < 1) {
        throw std::runtime_error(
            "cuMES binary: scientific fields are absent or invalid");
    }
    const std::size_t points = checked_product(
        static_cast<std::size_t>(equilibrium.ntheta),
        static_cast<std::size_t>(equilibrium.nzeta), "angular dimensions");
    const std::size_t half_size =
        checked_product(static_cast<std::size_t>(equilibrium.ns - 1), points,
                        "half-grid dimensions");
    const std::size_t full_size =
        checked_product(static_cast<std::size_t>(equilibrium.ns), points,
                        "full-grid dimensions");
    if (half_size > MAX_REAL_FIELD_ELEMENTS ||
        full_size > MAX_REAL_FIELD_ELEMENTS) {
        throw std::runtime_error(
            "cuMES binary: scientific fields exceed the resource cap");
    }
    for (auto& field : equilibrium.half_fields) {
        field = reader.read_f64_array(half_size, "half-grid field");
    }
    for (auto& field : equilibrium.full_fields) {
        field = reader.read_f64_array(full_size, "full-grid field");
    }

    skip_output_provenance(reader);
    read_input_params(reader, equilibrium);
    return equilibrium;
}

}  // namespace magnetic_coordinate
