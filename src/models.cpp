#include "../include/models.hpp"
#include "../include/surface.h"

#include <tuple>
#include <random>

/// Динамическая вязкость
double Air::dyn_vis(const double& T) const{
    return dyn_vis_c[0] + T * ( dyn_vis_c[1] + T * ( dyn_vis_c[2] + T * dyn_vis_c[3] ) );
}

/// Теплопроводность
double Air::lambda(const double& T) const{
    if(T > 10.){
        return lambda_c[0] * std::pow( T, lambda_c[1] ) / ( 1. + lambda_c[2] / T );
    }
    return lambda_c[0] * std::pow( 10., lambda_c[1] ) / ( 1. + lambda_c[2] / 10. );
}

/// Теплоемоксть молекулярного кислорода
double Air::cp_o2(const double& T) const{
    double t1;
    if(T < 200.){
        t1 = 1/(1.+o2_gas_spec_heat[1]*200.);
    }else if(T > 3000.){
        t1 = 1/(1.+o2_gas_spec_heat[1]*3000.);
    }else{
        t1 = 1/(1.+o2_gas_spec_heat[1]*T);
    }
    return o2_gas_spec_heat[0] + o2_gas_spec_heat[2]*t1*t1 + o2_gas_spec_heat[3]*t1*t1*t1 + o2_gas_spec_heat[4]*t1*t1*t1*t1;
}

/// Теплоемоксть молекулярного азота
double Air::cp_n2(const double& T) const{
    double t1;
    if(T < 200.){
        t1 = 1/(1.+n2_gas_spec_heat[1]*200.);
    }else if(T > 3000.){
        t1 = 1/(1.+n2_gas_spec_heat[1]*3000.);
    }else{
        t1 = 1/(1.+n2_gas_spec_heat[1]*T);
    }
    return n2_gas_spec_heat[0] + n2_gas_spec_heat[2] *t1*t1 + n2_gas_spec_heat[3] *t1*t1*t1 + n2_gas_spec_heat[4] *t1*t1*t1*t1;
}

double Dukowicz::get_cp_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p){
    return f*p.props->cp_vap[T] + o*env->cp_o2(T) + a*env->cp_n2(T) + ppp*env->cp_n2(T);
}

double Dukowicz::get_lambda_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p){
    fmol  = env->M / p.props->M * f;
    omol  = env->M / (32.0*1e-3) * o;
    aimol = env->M / (28.0*1e-3) * a;
    pmol = env->M / (30.0*1e-3) * ppp;

    return fmol*p.props->lambda_vap[T] + omol*env->lambda(T) + aimol*env->lambda(T) + pmol*env->lambda(T);
}

double Dukowicz::get_dyn_vis_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p){
    fmol = env->M/p.props->M * f;
    omol  = env->M/(32.0*1e-3) *o;  // o - конц. кислорода = 0.23
    aimol = env->M/(28.0*1e-3) *a;  // a - конц. азота = 0.77
    pmol = env->M/(30.0*1e-3) *ppp;

    double viso = env->dyn_vis(T);
    double visi = viso;
    double visp = viso;

    double visf = p.props->dyn_vis_vap[T];
    
    return fmol*visf+omol*viso+aimol*visi+pmol*visp;
}

void Dukowicz::solve(Particle& p,  Cell& c){
    T = (p.TGet() + c.TGet())/2.;

    f_vs =  p.props->p_sat[p.TGet()] / c.pGet();
    Yvs = p.props->M / ( p.props->M  + (env->M * ( ( 1/f_vs ) - 1 ) ) );
    By = (Yvs - c.YvGet()) / (1 - Yvs);

    o = 0.23 * (1.- 0.5 *Yvs);
    a = 0.77 * (1.- 0.5 * Yvs);

    dyn_vis_film = get_dyn_vis_film(T,0.5*Yvs,o,a,0,p);
    lambda_film = get_lambda_film(T,0.5*Yvs,o,a,0,p);
    cp_film = get_cp_film(T,0.5*Yvs,o,a,0,p);

    Pr = dyn_vis_film * cp_film / lambda_film;
    Re = ( c.vGet() +p.v_pulseGet() - p.vGet() ).length() * p.dGet()  * c.roGet() / dyn_vis_film;

    Nu = 2 + 0.6 * std::sqrt(Re) * std::pow(Pr, 1./3);
    q  = p.dGet() * pi * lambda_film * Nu * ( c.TGet() - p.TGet() );

    o = 0.23*(1.-Yvs);
    a = 0.77*(1.-Yvs);

    f_q = -By / 
        (
            ( get_cp_film(c.TGet(),0,0.23,0.77,0,p)  * c.TGet() ) - 
            ( get_cp_film(p.TGet(),Yvs,o,a,0,p) * p.TGet() ) -      
            (
                ( ( p.props->cp_vap[p.TGet()] * p.TGet() ) - ( get_cp_film(p.TGet(),0,0.23,0.77,0,p) * p.TGet() ) ) * 
                ( c.YvGet() - Yvs )
            )
        ) * C1;

    dT = (q * ( 1 + ( p.props->L[p.TGet()] * f_q ) ) / p.props->cp[p.TGet()] / p.mGet()) * p.dtGet();
    dm = (q * f_q) * p.dtGet();

    p.evaporation(dT, dm);

    // Передача источников тепла и массы
    c.recvHeat(q*p.nGet());
    c.recvMass(dm*p.nGet());
};

void Wave::solve(Particle& p,  Cell& c){
    // Веббер жидкостный
    We1 = p.props->ro[p.TGet()] * p.dGet() * std::pow( (c.vGet() +p.v_pulseGet() - p.vGet()).length(), 2) / p.props->sigma[p.TGet()] / 2;
    // Веббер газовый
    We2 = c.roGet() * p.dGet() * std::pow( (c.vGet() +p.v_pulseGet() - p.vGet()).length(), 2) / p.props->sigma[p.TGet()] / 2;
    Re = (c.vGet() +p.v_pulseGet() - p.vGet()).length() * p.dGet()  * p.props->ro[p.TGet()] / p.props->dyn_vis[p.TGet()] / 2;
    Oh = std::sqrt(We1) / Re;
    T_ = Oh * std::sqrt(We2);
    
    Lambda = 9.02 * (p.dGet() / 2) * (1 + 0.45*std::sqrt(Oh)) * (1 + 0.4 * std::pow(T_, 0.7)) / 
        std::pow((1 + 0.87 * std::pow(We2, 1.67)), 0.6
    );
    
    Omega = (0.34 + 0.38 * std::pow(We2, 1.5)) / 
        (1 + Oh) / 
        (1 + 1.4 * std::pow(T_, 0.6)) / 
        std::sqrt( p.props->ro[p.TGet()] * std::pow((p.dGet() / 2), 3) / 
        p.props->sigma[p.TGet()] 
    );
    tau_a = 3.726 * C2 * (p.dGet() / 2) / Lambda / Omega;

    dd = (p.dGet() - 2 * C1 * Lambda) / tau_a * p.dtGet();
    
    p.breakUp(dd);
}

Vec Momentum::drag(Particle& p,  Cell& c){
    Re = (c.vGet() +p.v_pulseGet() - p.vGet()).length() * p.dGet()  * c.roGet() / env->dyn_vis(c.TGet());
    if(Re < 1e3){
        Cd = 24 * (1 + 0.15*std::pow(Re, 0.687)) / Re;
    }else{
        Cd = 0.44;
    }
    return (c.vGet() +p.v_pulseGet() - p.vGet()) * (c.vGet() +p.v_pulseGet() - p.vGet()).length() * 3 * Cd * c.roGet() / (4 * p.props->ro[p.TGet()] * p.dGet());
}

std::tuple<Surface*, Vec> Momentum::checkCross(Particle& p, std::vector<Surface>& ss_vec, int& iskip){
    // Перевод координат из Vec в std::vector<double>(3)
    x0 = p.cGet().to_std(); // Начальная точка
    x1v = p.cGet() + (p.vGet() * p.dtGet()); // Шаг
    x1 = (x1v).to_std(); // Конечная точка
    
    int ielem = -1;

    for(Surface& ss: ss_vec){
        std::tie(xn,ielem) = ss.check_cross(x0,x1,iskip); // Проверка пересечения
        if (ielem != -1) {
            iskip = ielem;
            xnv = {xn[0],xn[1],xn[2]};
            return {&ss, xnv}; // Пересечение есть - {указатель на поверхность, координата}
        }
    }
    return {nullptr, {0,0,0}}; // Пересечения нет
};

void Momentum::wallInteraction(Particle& p, Surface& ss, int& ielem){
    std::vector<double> nn = ss.get_n(ielem);
    std::vector<double> vel = p.vGet().to_std();
    double cosa = ss.veccosa(vel,nn); // angle between vector and surface normal
    double sina = std::sqrt(1.0 - cosa*cosa);
    double vn = ss.vecmod(vel)*cosa; // module of normal velocity

    std::vector<double> vref = ss.reflect_ray_full(vel,nn);  // calculate reflected velocity vector
    double Wbin = p.props->ro[p.TGet()]*p.dGet()*vn*vn/p.props->sigma[p.TGet()]; // calculate in We number based on normal velocity
    double Wbout;

    std::vector<double> vres(3,0); // end velocity vector
    double refn; // coefficient for normal velocity
    
    // update velocity
    if (Wbin < 80) {  // Wbin < 80
        Wbout = 0.687*Wbin*std::exp(-0.04415*Wbin); // calculate reflected We number
        refn = std::sqrt(Wbout/Wbin);  // calculate coefficient for normal velocity
        for (int j = 0; j < 3; j++) {
            vres[j] = vref[j+3]*refn + vref[j+6]; // define reflection velocity based on normal and tang velocities
        }
    } else {
        double k = sin_alpha_k_table.k(sina); // take k from table
        std::vector<double> vt(3,0);
        for (int j = 0; j < 3; j++) { vt[j] = vref[j+6]; }
        double rnd = dis0(gen); // generate rnd number between 0..1
        double phirot = -(pi/k)*std::log(1.0-rnd*(1.0-std::exp(-k))); // calculate rotation angle for Vt based in rnd number
        vt = ss.vecrot_surf(vt,nn,phirot);  // rotate tang velocity on angle phirot on surface nn
        for (int j = 0; j < 3; j++) {
            vres[j] = vref[j+3] + vt[j]; // define reflection velocity based on normal and tang velocities
        }
    }

    // update droplet diameter
    if (Wbin >= 50 && Wbin < 300) {
        double dd = p.dGet()*(1 - 0.416*std::pow(10,-0.00102*Wbin));
        p.breakUp(dd);
    } else if (Wbin > 300) {
        double dd = p.dGet() * (1 - 0.2);
        p.breakUp(dd);
    }

    // update coordinate and velocity
    p.vSet() = {vres[0],vres[1],vres[2]};
};

void TurbulentDispersion::solve(Particle& p, Cell& c){
    // Модель расчитывает турбулентную скорость для частицы и замораживает ее на время жизни вихря.
    // По ситечении времени производится пересчет.
    if(p.time_pulseGet() <= 0){
        u = (p.vGet()).length()* 0.05;
        k = 1.5 * u * u;
        eps = 1e4 * k;

        u_ = std::sqrt(2.*k/3.);
        p.v_pulseSet() = {u_ * normal_dist(gen), u_ * normal_dist(gen), u_ * normal_dist(gen)};
        p.time_pulseSet() = std::min(k/eps, Ct*std::pow(k,1.5)/eps/(std::abs(c.vGet().length() + u_ - p.cGet().length())));
    }
};