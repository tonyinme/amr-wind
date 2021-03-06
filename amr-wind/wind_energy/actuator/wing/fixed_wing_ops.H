#ifndef FIXED_WING_OPS_H
#define FIXED_WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/FixedWing.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace ops {

template <>
struct ReadInputsOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.get("epsilon", wdata.eps_inp);
        pp.get("pitch", wdata.pitch);

        pp.get("airfoil_table", wdata.airfoil_file);
        pp.query("airfoil_type", wdata.airfoil_type);
        pp.queryarr("span_locs", wdata.span_locs);
        pp.queryarr("chord", wdata.chord_inp);
        AMREX_ALWAYS_ASSERT(wdata.span_locs.size() == wdata.chord_inp.size());

        amrex::Real max_chord =
            *std::max_element(wdata.chord_inp.begin(), wdata.chord_inp.end());
        amrex::Real max_eps =
            *std::max_element(wdata.eps_inp.begin(), wdata.eps_inp.end());
        amrex::Real search_radius = max_eps * max_chord * 3.0;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // clang-format off
        info.bound_box = amrex::RealBox(
            amrex::min(p1.x(), p2.x()) - search_radius,
            amrex::min(p1.y(), p2.y()) - search_radius,
            amrex::min(p1.z(), p2.z()) - search_radius,
            amrex::max(p1.x(), p2.x()) + search_radius,
            amrex::max(p1.y(), p2.y()) + search_radius,
            amrex::max(p1.z(), p2.z()) + search_radius
        );
        // clang-format on
    }
};

template <>
struct InitDataOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data)
    {
        wing::init_data_structures(data.meta(), data.grid());

        auto& meta = data.meta();
        {
            const int npts = data.meta().num_pts;
            const auto& grid = data.grid();
            // Process chord information
            const auto wlen = vs::mag(grid.pos.back() - grid.pos.front());
            RealList wx(npts);
            for (int i = 0; i < npts; ++i) {
                wx[i] = vs::mag(grid.pos[i]) / wlen;
            }
            meta.chord.resize(npts);
            ::amr_wind::interp::linear_monotonic(
                meta.span_locs, meta.chord_inp, wx, meta.chord);
        }

        meta.aflookup =
            AirfoilLoader::load_airfoil(meta.airfoil_file, meta.airfoil_type);
    }
};

template <>
inline const AirfoilTable& airfoil_lookup<FixedWing>(FixedWing::DataType& data)
{
    return *data.meta().aflookup;
}

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* FIXED_WING_OPS_H */
