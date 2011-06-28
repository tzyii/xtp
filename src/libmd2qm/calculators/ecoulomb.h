#ifndef _CALC_ESTATICS_H
#define	_CALC_ESTATICS_H

#include "qmpair.h"
#include "qmcalculator.h"

/**
    \brief Coulomb contribution to site energies

Callname: ecoulomb

Electrostatic contributions to site-energies of \f$E_i\f$ are computed from Coulomb's law for all sites \f$i\f$
    \f[E_{i}  = \frac{1}{4 \pi \epsilon_0} \sum_{a_i} \sum_{b_k, k\neq i }
\frac{ \left( q^c_{a_i} - q^n_{a_i} \right) q^n_{b_k}}{ \epsilon_{s} r_{a_i b_k}}\f]
where \f$r_{a_i b_k}=|\vec{r}_{a_i} - \vec{r}_{b_k}|\f$
is the distance between atoms \f$a_i\f$ (which are in site \f$i\f$) and \f$b_k\f$ (which are in site \f$k\f$). The minimum image convention is used and the partial charges are taken for neutral, \f$q^n\f$, or charged, \f$q^c\f$, sites. Finally \f$\epsilon_s\f$ is the static relative dielectric constant. Depending on the option chosen in the main file \f$\epsilon_s\f$ is either constant or distance-dependent. In the latter case a screening length must be specified \f$s\f$.

*/
class CalcEstatics : public QMCalculator
{
public:
    CalcEstatics() {};
    ~CalcEstatics() {};

    const char *Description() { return "Electrostatic contribution to site-energies"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

    double CalcPot(Topology *atop, Molecule *mol);
    /// distance dependent dielectric constant
    double dist_dep_eps(const double &dist);
    /// constant epsilon
    double constant_epsilon(const double &dist);

private:
    // dielectric constant
    double _epsilon_dielectric;
    // parameter describing the decay of eps
    double _s_eps;
    Property * _options;
    // function pointer 
    double (CalcEstatics::*_estatic_method)(const double &);
};

#endif	/* _CALC_ESTATICS_H */
