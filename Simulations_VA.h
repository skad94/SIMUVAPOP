#pragma once
//#include "Simulations_VA.h"
#define PI_nb 3.14159

class continuousProcessSimulable
{
private:
    double m_drift;
    double m_diffusion;
    double m_start;

public:
    continuousProcessSimulable();
    continuousProcessSimulable(double dri, double diffusion, double start)
        :m_drift(dri), m_diffusion(diffusion), m_start(start)
    {}
    double getStart() const { return m_start; };
    double getDrift() const { return m_drift; };
    double getDiffusion() const { return m_diffusion; };
    std::vector<double> fillPath(const std::vector<double>& data, const size_t& finalSize) const;
    std::vector<double> Simul_Brownian(const int& sizeSample, double mu = 0.0, double sigma = 1.0) const;
    virtual std::vector<double> simulPath(const int& sizeSample) const = 0;
    // virtual std::vector<double> SimulPath(const int& sizeSample, const std::vector<double>& brownianPath) const = 0;
};

class PoissonProcess
{
private:
    double m_intensity;
public:
    PoissonProcess(double intensity);
    double getIntensity() const { return m_intensity; }
    virtual std::vector<double> simulPath(const int& eventAmount) const;
};

class BS_Process : public continuousProcessSimulable
{
public:
    BS_Process(double dri, double diffusion, double start)
        :continuousProcessSimulable(dri, diffusion, start)
    {}
    virtual std::vector<double> simulPath(const int& sizeSample) const;
    std::vector<double> BS_Simul_Euler(const int& sizeSample) const;
    std::vector<double> BS_Simul_Euler(const int& sizeSample, const std::vector<double>& brownianPath) const;
    std::vector<double> BS_Simul_milstein(const int& sizeSample) const;
    std::vector<double> BS_Simul_milstein(const int& sizeSample, const std::vector<double>& brownianPath) const;
};

class vasicek_Process : public continuousProcessSimulable
{
private:

    double m_speed_reversion;
    double m_mean_reversion;
public:
    vasicek_Process(double dri, double diffusion, double start, double speed_reversion, double mean_reversion)
        :continuousProcessSimulable(dri, diffusion, start), m_speed_reversion(speed_reversion), m_mean_reversion(mean_reversion)
    {}
    double getSpeed_reversion() const { return m_speed_reversion; }
    double getMean_reversion() const { return m_mean_reversion; }
    virtual std::vector<double> simulPath(const int& sizeSample) const;
    std::vector<double> vasicek_Simul_Euler(const int& sizeSample) const;
    //std::vector<double> vasicek_Simul_Euler(const int& sizeSample, const std::vector<double>& brownianPath) const;
    std::vector<double> vasicek_Simul_milstein(const int& sizeSample) const;
    //std::vector<double> vasicek_Simul_milstein(const int& sizeSample, const std::vector<double>& brownianPath) const;
};
std::vector<double> Simul_Brownien(int nb, double mu, double sigma);
std::vector<double> Markov_chain_3(int nb, int depart, double p11, double p12, double p22, double p23, double p33, double p31);
void Markov_chain(std::vector<std::vector<double>> p, int nb_etat, int taille, int depart);


std::vector<double> Proc_CIR_Simul(int nb, double origine, double moyenne, double diffusion, double speed_retour);
std::vector<double> Proc_pseudo_CEV(int nb, double origine, double moyenne, double nu, double puissance);
double sd_hat(std::vector<double> X, int n);
std::vector<double> Trier(std::vector<double> tab, int taille);
int TVI(std::vector<double> t, int taille, double e);

std::vector<double> BS_Simul(int nb, double origine, double mu, double diffusion);

std::vector<double> BS_Simul_milstein(int nb, double origine, double mu, double diffusion);
double f(double x, double lambda);
double h(double x);

void Simul_Sample(std::vector<double> val, int size_sample, double sd, double su, double p);
void display(const std::vector<double>& vect);
double MAX(const double& A, const double& B);

double Frequence(const std::vector<double>& vect, double s);


//void Gain_cumule(std::vector<double> x, int nb, double x0, double sa, double sb, double p);
//void M_C_Gain_cumule(std::vector<double> x, int nb, int nb_mc, double x0, double sa, double sb, double p);
//void M_C_Gain_cumule_Max(std::vector<double> x, int nb, int nb_mc, int maxx, double x0, double sa, double sb, double p);
//void Fcts_Barriere_up_and_down(std::vector<double> val, int taille, double sup, double inf);
//std::vector<double> Temps_Atteinte(int nb, int nb_stop, double x0, double sa, double sb, double p);
void Call_BS_MC(int nb_MC, int nb_pas, double spot, double mu, double diffusion, double K);