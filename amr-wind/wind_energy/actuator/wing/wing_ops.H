#ifndef WING_OPS_H
#define WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/core/vs/vector_space.H"

#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_AmrCore.H"


namespace amr_wind {
namespace actuator {
namespace utils {
class ActParser;
}

namespace wing {

/** Parse common data for modeling fixed wings as actuators
 */
void read_inputs(WingBaseData&, ActInfo&, const utils::ActParser&);

/** Initialize core data structures when modeling fixed wings
 */
void init_data_structures(WingBaseData&, ActGrid&);

void prepare_netcdf_file(
    const std::string&, const WingBaseData&, const ActInfo&, const ActGrid&);

void write_netcdf(
    const std::string&,
    const WingBaseData&,
    const ActInfo&,
    const ActGrid&,
    const amrex::Real);

} // namespace wing

namespace ops {

template <typename T>
const typename AirfoilTraits<T>::AirfoilLookup&
airfoil_lookup(typename T::DataType& data);

template <typename ActTrait, typename SrcTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_pos>::type>
{
    void operator()(typename ActTrait::DataType&) {}
};

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_vel>::type>
{
    void operator()(typename ActTrait::DataType& data) {
    // Check if FLLC is active
    auto& meta = data.meta(); 
    if (! meta.fllc_use || ! meta.fllc_is_init) return;
    // Implement the FLLC
    FLLCOp<ActTrait, ActSrcLine> op(data);
    op(data);
    }
};

template <typename ActTrait>
struct ComputeForceOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_force>::type>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        const auto& chord = wdata.chord;
        const auto& aflookup = airfoil_lookup<ActTrait>(data);

        amrex::Real total_lift = 0.0;
        amrex::Real total_drag = 0.0;
        for (int ip = 0; ip < npts; ++ip) {
            const auto& tmat = grid.orientation[ip];
            // Effective velocity at the wing control point in local frame
            auto wvel = tmat & grid.vel[ip];

            // Set spanwise component to zero to get a pure 2D velocity
            wvel.y() = 0.0;

            const auto vmag = vs::mag(wvel);
            const auto aoa = std::atan2(wvel.z(), wvel.x());

            // Make up some Cl, Cd values
            amrex::Real cl, cd;
            aflookup(aoa, cl, cd);

            // Assume unit chord
            const auto qval = 0.5 * vmag * vmag * chord[ip] * dx[ip];
            const auto lift = qval * cl;
            const auto drag = qval * cd;
            // Determine unit vector parallel and perpendicular to velocity
            // vector
            const auto drag_dir = wvel.unit() & tmat;
            const auto lift_dir = drag_dir ^ tmat.y();

            // Compute force on fluid from this section of wing
            grid.force[ip] = -(lift_dir * lift + drag * drag_dir);

            // Assign values for output
            wdata.vel_rel[ip] = wvel;
            wdata.aoa[ip] = amr_wind::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;

        // Initialize the fllc flag
        wdata.fllc_is_init = true;

    }
};

template <typename ActTrait, typename SrcTrait>
struct ProcessOutputsOp<
    ActTrait,
    SrcTrait,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::process_outputs>::type>
{
private:
    typename ActTrait::DataType& m_data;

    std::string m_out_dir{""};
    std::string m_nc_filename{""};
    int m_out_freq{10};

public:
    ProcessOutputsOp(typename ActTrait::DataType& data) : m_data(data) {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        wing::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        wing::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

template <typename ActTrait, typename SrcTrait>
struct FLLCOp<
    ActTrait, 
    SrcTrait,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_vel>::type>
{

    //! The velocity correction from filtered lifting line theory
    VecList u_les_fllc;
    VecList u_opt_fllc;
    VecList G_fllc;
    VecList dG_fllc;

    FLLCOp(typename ActTrait::DataType& data){

        auto& wdata = data.meta();
        const int npts = wdata.num_pts;

         G_fllc.assign(npts, vs::Vector::zero());
         dG_fllc.assign(npts, vs::Vector::zero());
         u_les_fllc.assign(npts, vs::Vector::zero());
         u_opt_fllc.assign(npts, vs::Vector::zero());
        
    }

    void operator()(typename ActTrait::DataType& data) {
        
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        auto& du_fllc = wdata.du_fllc;

        // Step 1
        // Compute G
        for (int ip = 0; ip < npts; ++ip) {

            const auto& tmat = grid.orientation[ip];
            auto force = grid.force[ip];
//            auto vel =  grid.vel[ip];
            auto vel =  wdata.vel_rel[ip] & tmat;
            auto dx =  wdata.dx[ip];
            const auto vmag = vs::mag(vel);
            const auto vmag2 = vmag * vmag;

            const auto fv = force & vel;

            G_fllc[ip] = (force - vel * fv / vmag2) / dx;

            }

        // Step 2
        // Compute the gradient of the lift dG
       dG_fllc[0] = G_fllc[0];
       dG_fllc[npts-1] = -1 * G_fllc[npts-1];
       for (int ip = 1; ip < npts-1; ++ip) {

            dG_fllc[ip] = 0.5 * (G_fllc[ip+1] - G_fllc[ip - 1]);

          }

        // Step 3
        // Copmute the induced velocity from epsilon LES and
        // from the optimal epsilon
        for (int ip = 0; ip < npts; ++ip) {

            const auto eps_les = wdata.eps_inp[0];
            const auto eps_opt = wdata.epsilon_chord[0] * wdata.chord[ip];
                    
            for (int jp = 0; jp < npts; ++jp) {

                if (ip==jp) continue; 

                const auto dr = vs::mag(grid.vel_pos[ip] - grid.vel_pos[jp]);
                const auto& vel =  wdata.vel_rel[jp];
                const auto vmag = vs::mag(vel);

                auto coefficient = 1.0 / (-4.0 * amr_wind::utils::pi() * dr * vmag);
                const auto cLes = 1.0 - std::exp(-dr * dr / (eps_les * eps_les));
                const auto cOpt = 1.0 - std::exp(-dr * dr / (eps_opt * eps_opt));

                // The sign of the induced velocity depends on which side of the blade we are on
                if (ip<jp) coefficient*=-1;
           
               u_les_fllc[ip] = u_les_fllc[ip] - dG_fllc[jp] * coefficient * cLes;
               u_opt_fllc[ip] = u_opt_fllc[ip] - dG_fllc[jp] * coefficient * cOpt;
                
            }

            // Relaxation to compute the induced velocity
            const double f = 0.1;
            du_fllc[ip] = (1.-f) * du_fllc[ip] + f * (u_opt_fllc[ip] - u_les_fllc[ip]);
            
        }

        // Step 4
        // Update the velocity using the previous time
        for (int ip = 0; ip < npts; ++ip) {

            // Update the velocity with the correction from the previous time
            grid.vel[ip] = grid.vel[ip] + du_fllc[ip];
        }

    }

};

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* WING_OPS_H */
