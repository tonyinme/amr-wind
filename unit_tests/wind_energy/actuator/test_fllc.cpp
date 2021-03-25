#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"


namespace amr_wind_tests {
namespace {
class ActFLLCtest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{32.0, 32.0, 32.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace act = amr_wind::actuator;
namespace vs = amr_wind::vs;



} // namespace

} // namespace amr_wind_tests
