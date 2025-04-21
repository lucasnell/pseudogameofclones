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
    WaspPop wasps;                  // wasps / mummies in this field
    bool empty;                     // whether no aphids are in this field
    double pred_rate;               // predation on aphids and mummies
    double K;                       // unparasitized aphid carrying capacity
    double K_y;                     // parasitized aphid carrying capacity
    double z = 0;                   // sum of all living aphids at time t
    double S = 0;                   // effect of density dependence on aphids
    double S_y = 0;                 // effect of dd on parasitized aphids
    uint32 n_fields;                // total # fields
    double extinct_N;               // threshold for calling something extinct
    bool constant_wasps;            // keep wasp abundance constant?



    OneField()
        : aphids(), wasps(), empty(true), pred_rate(0),
          K(0), K_y(1), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types
     (alate vs apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are
     alate/apterous/parasitized.
     */
    OneField(const std::vector<AphidPop>& aphids_,
             const WaspPop& wasps_,
             const double& K_,
             const double& K_y_,
             const double& pred_rate_,
             const double& extinct_N_,
             const bool& constant_wasps_)
        : aphids(aphids_),
          wasps(wasps_),
          empty(true),
          pred_rate(pred_rate_),
          K(K_),
          K_y(K_y_),
          extinct_N(extinct_N_),
          constant_wasps(constant_wasps_) {

        for (AphidPop& aphid : aphids) {
            double N = aphid.total_aphids();
            if (N < extinct_N) {
                aphid.clear();
            } else if (empty) empty = false;
        }

    };



    OneField(const OneField& other)
        : aphids(other.aphids), wasps(other.wasps),
          empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), K_y(other.K_y), z(other.z),
          S(other.S), S_y(other.S_y),
          extinct_N(other.extinct_N),
          constant_wasps(other.constant_wasps){};

    OneField& operator=(const OneField& other) {
        aphids = other.aphids;
        wasps = other.wasps;
        empty = other.empty;
        pred_rate = other.pred_rate;
        K = other.K;
        K_y = other.K_y;
        z = other.z;
        S = other.S;
        S_y = other.S_y;
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

    using iterator = std::vector<AphidPop>::iterator;
    using const_iterator = std::vector<AphidPop>::const_iterator;
    iterator begin() {
        return aphids.begin();
    }
    const_iterator begin() const {
        return aphids.begin();
    }
    iterator end() {
        return aphids.end();
    }
    const_iterator end() const {
        return aphids.end();
    }


    /*
     Clear to no aphids or mummies
     */
    void clear() {
        for (AphidPop& ap : aphids) ap.clear();
        wasps.mummies.clear();
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
        wasps.mummies.clear(surv);
        if (arma::accu(wasps.mummies.Y) < extinct_N) wasps.mummies.Y.fill(0);
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
        double tm = arma::accu(wasps.mummies.Y);
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
              const std::vector<AphidPop>& aphids_,
              const WaspPop& wasps_,
              const std::vector<double>& K,
              const std::vector<double>& K_y,
              const std::vector<double>& pred_rate,
              const double& extinct_N_,
              const std::vector<bool>& constant_wasps,
              const double& alate_field_disp_p_,
              const double& wasp_disp_m0_,
              const double& wasp_disp_m1_,
              const std::vector<double>& wasp_field_attract_)
        : eng(), fields(),
          alate_field_disp_p(alate_field_disp_p_),
          wasp_disp_m0(wasp_disp_m0_),
          wasp_disp_m1(wasp_disp_m1_),
          wasp_field_attract(wasp_field_attract_),
          extinct_N(extinct_N_),
          total_stages(0) {

        // This vector can't have values < 0:
        double wfa_min = *min_element(wasp_field_attract.begin(),
                                      wasp_field_attract.end());
        if (wfa_min < 0) stop("\nmin(wasp_field_attract) < 0.\n");
        // It can have zeros, but can't have all zeros:
        double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                         wasp_field_attract.end(), 0.0);
        if (wfa_sum <= 0) stop("\nwasp_field_attract sums to <= 0.\n");
        // This makes it sum to 1:
        for (double& x : wasp_field_attract) x /= wfa_sum;


        fields.reserve(n_fields);
        for (uint32 i = 0; i < n_fields; i++) {
            fields.push_back(
                OneField(aphids_, wasps_, K[i], K_y[i], pred_rate[i],
                         extinct_N, constant_wasps[i]));
            total_stages += 1;
            const OneField& field(fields[i]);
            total_stages += field.wasps.mummies.Y.n_elem;
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

    using iterator = std::vector<OneField>::iterator;
    using const_iterator = std::vector<OneField>::const_iterator;
    iterator begin() {
        return fields.begin();
    }
    const_iterator begin() const {
        return fields.begin();
    }
    iterator end() {
        return fields.end();
    }
    const_iterator end() const {
        return fields.end();
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
                double z, p_out, lz;
                for (OneField& field : fields) {
                    if (field.wasps.Y <= 0) continue;
                    z = field.total_aphids();
                    if (z > 0) {
                        lz = std::log(z);
                        p_out = wasp_disp_m0 * std::exp(-wasp_disp_m1 * lz);
                        if (p_out > 1) p_out = 1; // can happen if z < 1
                        from_wasp_pool += (field.wasps.Y * p_out);
                        field.wasps.Y *= (1 - p_out);
                    } else {
                        from_wasp_pool += field.wasps.Y;
                        field.wasps.Y = 0;
                    }
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
