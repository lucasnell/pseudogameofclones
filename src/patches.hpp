# ifndef __PSEUDOGAMEOFCLONES_PATCHES_H
# define __PSEUDOGAMEOFCLONES_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <cmath>                // std::exp, std::log
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "math.hpp"             // distributions
#include "aphids.hpp"           // aphid classes
#include "wasps.hpp"            // wasp classes



using namespace Rcpp;

// Necessary here to declare friendship
class AllFields;



// Info for a single perturbation for a field
struct PerturbInfo {

    // When to perturb.
    uint32 time;
    // Where to perturb.
    uint32 field;
    // Number to multiply abundance by.
    double multiplier;
    /*
     Who to perturb.
     Indices less than the # aphids points to the particular
     aphid line, equal to the # aphids is mummies, and greater than
     the # aphids is adult wasps.
     */
    uint32 index;

    PerturbInfo() : time(), field(), multiplier(), index() {}
    PerturbInfo(const uint32& when, const uint32& where,
                const uint32& who, const double& how)
        : time(when), field(where), multiplier(how), index(who) {}
    PerturbInfo& operator=(const PerturbInfo& other) {
        time = other.time;
        field = other.field;
        multiplier = other.multiplier;
        index = other.index;
        return *this;
    }

};











/*
 One field.

 For the constructors below, all vector arguments should have a length equal
 to the number of aphid lines.
*/
class OneField {

    friend class AllFields;

    /*
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i) {
        AphidPop& ap(aphids[i]);
        double N = ap.total_aphids();
        if (N < extinct_N) {
            ap.clear();
        } else {
            empty = false;
            ap.extinct = false;
        }
        return;
    }




public:

    std::vector<AphidPop> aphids;   // aphids in this field
    MummyPop mummies;               // mummies in this field
    WaspPop wasps;                  // wasps in this field
    bool empty;                     // whether no aphids are in this field
    double pred_rate;               // predation on aphids and mummies
    double K;                       // unparasitized aphid carrying capacity
    double K_y;                     // parasitized aphid carrying capacity
    double z = 0;                   // sum of all living aphids at time t
    double S = 0;                   // effect of density dependence on aphids
    double S_y = 0;                 // effect of dd on parasitized aphids
    uint32 n_fields;                // total # fields
    uint32 this_j;                  // index for this field
    double extinct_N;               // threshold for calling something extinct
    bool constant_wasps;            // keep wasp abundance constant?



    OneField()
        : aphids(), mummies(), wasps(), empty(true), pred_rate(0),
          K(0), K_y(1), this_j(0), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types
     (alate vs apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are
     alate/apterous/parasitized.
     */
    OneField(const double& sigma_x,
             const double& sigma_y,
             const double& rho,
             const bool& aphid_demog_error,
             const bool& wasp_demog_error,
             const arma::mat& attack_surv_,
             const double& K_,
             const double& K_y_,
             const std::vector<std::string>& aphid_name,
             const std::vector<arma::cube>& leslie_mat,
             const arma::cube& aphid_density_0,
             const std::vector<double>& alate_b0,
             const std::vector<double>& alate_b1,
             const std::vector<uint32>& field_disp_start,
             const std::vector<uint32>& living_days,
             const double& pred_rate_,
             const uint32& this_j_,
             const double& extinct_N_,
             const arma::vec& mum_density_0,
             const double& mum_smooth,
             const arma::vec& rel_attack_,
             const double& a_,
             const double& k_,
             const double& h_,
             const double& wasp_density_0_,
             const uint32& wasp_delay_,
             const double& sex_ratio_,
             const double& s_y_,
             const bool& constant_wasps_)
        : aphids(),
          mummies(mum_density_0, mum_smooth),
          wasps(rel_attack_, a_, k_, h_,
                wasp_density_0_, wasp_delay_,
                sex_ratio_, s_y_, sigma_y, wasp_demog_error),
          empty(true),
          pred_rate(pred_rate_),
          K(K_),
          K_y(K_y_),
          this_j(this_j_),
          extinct_N(extinct_N_),
          constant_wasps(constant_wasps_) {

        uint32 n_lines = aphid_name.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name[i], sigma_x, rho, aphid_demog_error,
                        attack_surv_.col(i),
                        leslie_mat[i], aphid_density_0.slice(i),
                        alate_b0[i], alate_b1[i], field_disp_start[i],
                        living_days[i]);
            aphids.push_back(ap);
            double N = aphids.back().total_aphids();
            if (N < extinct_N) {
                aphids.back().clear();
            } else if (empty) empty = false;
        }

    };



    OneField(const OneField& other)
        : aphids(other.aphids), mummies(other.mummies), wasps(other.wasps),
          empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), K_y(other.K_y), z(other.z),
          S(other.S), S_y(other.S_y),
          this_j(other.this_j),
          extinct_N(other.extinct_N),
          constant_wasps(other.constant_wasps){};

    OneField& operator=(const OneField& other) {
        aphids = other.aphids;
        mummies = other.mummies;
        wasps = other.wasps;
        empty = other.empty;
        pred_rate = other.pred_rate;
        K = other.K;
        K_y = other.K_y;
        z = other.z;
        S = other.S;
        S_y = other.S_y;
        this_j = other.this_j;
        extinct_N = other.extinct_N;
        constant_wasps = other.constant_wasps;
        return *this;
    }

    inline uint32 size() const noexcept {
        return aphids.size();
    }

    AphidPop& operator[](const uint32& idx) {
        return aphids[idx];
    }
    const AphidPop& operator[](const uint32& idx) const {
        return aphids[idx];
    }




    /*
     Clear to no aphids or mummies
     */
    void clear() {
        for (AphidPop& ap : aphids) ap.clear();
        mummies.clear();
        empty = true;
        return;
    }
    /*
     Clear some of the aphids and mummies
     */
    void clear(const double& surv) {
        empty = true;
        for (AphidPop& ap : aphids) {
            ap.clear(surv);
            double N = ap.total_aphids();
            if (N < extinct_N) {
                ap.clear();
            } else {
                empty = false;
                ap.extinct = false;
            }
        }
        mummies.clear(surv);
        if (arma::accu(mummies.Y) < extinct_N) mummies.Y.fill(0);
        return;
    }

    // Total (living) aphids in field
    inline double total_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_aphids();
        }
        return ta;
    }
    // Total UNparasitized aphids in field
    inline double total_unpar_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.apterous.total_aphids() + ap.alates.total_aphids();
        }
        return ta;
    }
    // Total mummies in field
    inline double total_mummies() const {
        double tm = arma::accu(mummies.Y);
        return tm;
    }

    /*
     Adult and juvenile numbers for alates and apterous:
     */
    inline double total_adult_apterous() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_adult_apterous();
        }
        return ta;
    }
    inline double total_juven_apterous() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_juven_apterous();
        }
        return ta;
    }
    inline double total_adult_alates() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_adult_alates();
        }
        return ta;
    }
    inline double total_juven_alates() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_juven_alates();
        }
        return ta;
    }





    /*
     Iterate one time step, after calculating dispersal numbers
     */
    void update(pcg32& eng);

    arma::mat remove_field_dispersers(const double& disp_prop);
    void add_field_dispersers(const arma::mat& D);


};











// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================



struct AllStageInfo {

    std::vector<uint32> field;
    std::vector<std::string> line;
    std::vector<std::string> type;
    std::vector<uint32> stage;
    std::vector<double> N;

    void reserve(const uint32& total_stages) {
        field.reserve(total_stages);
        line.reserve(total_stages);
        type.reserve(total_stages);
        stage.reserve(total_stages);
        N.reserve(total_stages);
        return;
    }

    void push_back(const uint32& field_,
                   const std::string& line_,
                   const std::string& type_,
                   const uint32& stage_,
                   const double& N_) {
        field.push_back(field_);
        line.push_back(line_);
        type.push_back(type_);
        stage.push_back(stage_);
        N.push_back(N_);
        return;
    }

    DataFrame to_data_frame() const {
        DataFrame out = DataFrame::create(
            _["field"] = field,
            _["line"] = line,
            _["type"] = type,
            _["stage"] = stage,
            _["N"] = N);
        return out;
    }
};




class AllFields {

    pcg32 eng;   // RNG

public:

    std::vector<OneField> fields;

    double alate_field_disp_p;
    double wasp_disp_m0;
    double wasp_disp_m1;
    std::vector<double> wasp_field_attract;
    double extinct_N;

    // Total stages for all wasps, mummies, and aphids. Used for output.
    uint32 total_stages;

    AllFields()
        : eng(), fields(),
          alate_field_disp_p(), wasp_disp_m0(), wasp_disp_m1(),
          wasp_field_attract(),
          extinct_N(), total_stages() {};

    AllFields(const AllFields& other)
        : eng(other.eng),
          fields(other.fields),
          alate_field_disp_p(other.alate_field_disp_p),
          wasp_disp_m0(other.wasp_disp_m0),
          wasp_disp_m1(other.wasp_disp_m1),
          wasp_field_attract(other.wasp_field_attract),
          extinct_N(other.extinct_N),
          total_stages(other.total_stages){};

    AllFields& operator=(const AllFields& other) {
        eng = other.eng;
        fields = other.fields;
        alate_field_disp_p = other.alate_field_disp_p;
        wasp_disp_m0 = other.wasp_disp_m0;
        wasp_disp_m1 = other.wasp_disp_m1;
        wasp_field_attract = other.wasp_field_attract;
        extinct_N = other.extinct_N;
        total_stages = other.total_stages;
        return *this;
    };


    AllFields(const uint32& n_fields,
              const std::vector<uint32>& wasp_delay,
              const double& sigma_x,
              const double& sigma_y,
              const double& rho,
              const bool& aphid_demog_error,
              const bool& wasp_demog_error,
              const std::vector<double>& K,
              const std::vector<double>& K_y,
              const arma::mat& attack_surv,
              const std::vector<std::string>& aphid_name,
              const std::vector<arma::cube>& leslie_mat,
              const std::vector<arma::cube>& aphid_density_0,
              const std::vector<double>& alate_b0,
              const std::vector<double>& alate_b1,
              const std::vector<uint32>& field_disp_start,
              const std::vector<uint32>& living_days,
              const std::vector<double>& pred_rate,
              const double& extinct_N_,
              const arma::mat& mum_density_0,
              const double& mum_smooth,
              const arma::vec& rel_attack,
              const double& a,
              const double& k,
              const double& h,
              const std::vector<double>& wasp_density_0,
              const double& sex_ratio,
              const std::vector<double>& s_y,
              const std::vector<bool>& constant_wasps,
              const double& alate_field_disp_p_,
              const double& wasp_disp_m0_,
              const double& wasp_disp_m1_,
              const std::vector<double>& wasp_field_attract_,
              const std::vector<uint64>& seeds)
        : eng(), fields(),
          alate_field_disp_p(alate_field_disp_p_),
          wasp_disp_m0(wasp_disp_m0_),
          wasp_disp_m1(wasp_disp_m1_),
          wasp_field_attract(wasp_field_attract_),
          extinct_N(extinct_N_),
          total_stages(0) {

        seed_pcg(eng, seeds);

        //' Make sure `wasp_field_attract` sums to 1 (negative values
        //' and a sum <= 0 are already checked for in sim_pseudogameofclones_cpp):
        double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                         wasp_field_attract.end(), 0.0);
        for (double& x : wasp_field_attract) x /= wfa_sum;

        fields.reserve(n_fields);
        for (uint32 i = 0; i < n_fields; i++) {
            fields.push_back(
                OneField(sigma_x, sigma_y, rho,
                            aphid_demog_error, wasp_demog_error,
                         K[i], K_y[i],
                         attack_surv,
                         aphid_name, leslie_mat, aphid_density_0[i],
                         alate_b0, alate_b1,
                         field_disp_start,
                         living_days, pred_rate[i],
                         extinct_N, mum_density_0, mum_smooth,
                         rel_attack, a, k, h,
                         wasp_density_0[i], wasp_delay[i],
                         sex_ratio,
                         s_y[i], constant_wasps[i], eng));
            total_stages += 1;
            const OneField& field(fields[i]);
            total_stages += field.mummies.Y.n_elem;
            for (const AphidPop& aphid : field.aphids) {
                total_stages += aphid.apterous.X.n_elem;
                total_stages += aphid.alates.X.n_elem;
                total_stages += aphid.paras.X.n_elem;
            }
        }

    };

    inline uint32 size() const noexcept {
        return fields.size();
    }

    OneField& operator[](const uint32& idx) {
        return fields[idx];
    }
    const OneField& operator[](const uint32& idx) const {
        return fields[idx];
    }

    uint32 n_lines() const {
        return fields[0].aphids.size();
    }

    void reseed(const std::vector<uint64>& seeds) {
        seed_pcg(eng, seeds);
        return;
    }


    // Alate dispersal across fields:
    void across_field_disp_alates() {
        if (fields.size() > 1 && alate_field_disp_p > 0) {
            arma::mat D = fields[0].remove_field_dispersers(alate_field_disp_p);
            for (uint32 i = 1; i < fields.size(); i++) {
                D += fields[i].remove_field_dispersers(alate_field_disp_p);
            }
            D /= static_cast<double>(fields.size());
            for (OneField& field : fields) field.add_field_dispersers(D);
        }
        return;
    }

    // Wasp dispersal across fields:
    void across_field_disp_wasps() {
        if (fields.size() > 1 && wasp_disp_m0 > 0) {
            double from_wasp_pool = 0;
            if (wasp_disp_m1 != 0) {
                double p_out, lz;
                for (OneField& field : fields) {
                    lz = std::log(field.total_aphids());
                    p_out = wasp_disp_m0 * std::exp(-wasp_disp_m1 * lz);
                    from_wasp_pool += (field.wasps.Y * p_out);
                    field.wasps.Y *= (1 - p_out);
                }
            } else {
                for (OneField& field : fields) {
                    from_wasp_pool += (field.wasps.Y * wasp_disp_m0);
                    field.wasps.Y *= (1 - wasp_disp_m0);
                }
            }
            for (uint32 i = 0; i < fields.size(); i++) {
                OneField& field(fields[i]);
                field.wasps.Y += (from_wasp_pool * wasp_field_attract[i]);
            }
        }
        return;
    }

    void do_perturb(std::deque<PerturbInfo>& perturbs, const uint32& t);


    // Update for one time step.
    // Returns true if all fields are empty
    bool update(const uint32& t,
                std::deque<PerturbInfo>& perturbs);


    AllStageInfo out_all_info() const;

    // It's assumed this vector is in the same order as the outputs are above!
    // Make sure this happens from the R side.
    void from_vector(std::vector<double>& N);

    // Set new parameters
    void set_new_pars(const std::vector<double>& K_,
                      const std::vector<double>& alate_b0_,
                      const std::vector<double>& alate_b1_,
                      const double& alate_field_disp_p_,
                      const std::vector<double>& K_y_,
                      const std::vector<double>& s_y_,
                      const double& a_,
                      const double& k_,
                      const double& h_,
                      const double& wasp_disp_m0_,
                      const double& wasp_disp_m1_,
                      const std::vector<double>& wasp_field_attract_,
                      const double& mum_smooth_,
                      const std::vector<double>& pred_rate_);



};




#endif
