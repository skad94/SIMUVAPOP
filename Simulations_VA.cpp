#include <vector>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>
#include "Simulation.h"
//std::vector<double> PoissonProcess::simulPath(const int& y) const { return std::vector<double>(y); }
//PoissonProcess::PoissonProcess(double a) : ProcessSimulable(), m_intensity(a) {}
double sumSubArray(const std::vector<double>& data, const int& start, const int& stop)
{
    double res = data[start];
    for (int i = start + 1; i < stop; ++i)
        res += data[i];
    return res;
}
double 
unif()
{
    // return [0,1] uniform realisation
    return rand() / (RAND_MAX + 1.0);
}
ProcessSimulable::ProcessSimulable()
    :m_start(0.), m_drift(0.), m_diffusion(0.)
{}

std::vector<double> 
ProcessSimulable::fillPath(const std::vector<double>& data, const size_t& finalSize) const
{
    //in order to generate a finalSize's length vector from data.
    // linear interpolation between points in data
    std::vector<double> res(finalSize);
    int coarsityRatio = finalSize / data.size();
    for (int tk = 0; tk < finalSize; ++tk)
    {
        int tmp = tk % finalSize;
        int idxInData = tk / data.size();
        // interpolation between two times from data that surround the time we want to add
        res[tk] = (double)(tmp / coarsityRatio) * data[idxInData + 1] + ((double)(coarsityRatio - tmp) / (double)(coarsityRatio)) * data[idxInData];
    }
    return res;
}

std::vector<double>
ProcessSimulable::Simul_Brownian(const int& sizeSample, double mu, double sigma) const
{
    ///Return a Brownian path with possible drift (mu) and diffusion coefficient (sigma)
    std::vector<double> W(sizeSample,0.); // brownian
    /// Sample normal centered and reduced via box Muller method
    for (int i = 1; i < sizeSample; ++i)
    {
        double epsilon = sqrt(-2 * log(unif())) * cos(2 * PI_nb * unif());// Normal 
        epsilon = sigma * epsilon + mu;
        W[i] = W[i - 1] + sqrt(1.0 / (double)(sizeSample)) * epsilon;
    }
    return W;
}

//std::vector<double> 
//ProcessSimulable::simulPath(const int& sizeSample) const
//{
//    std::vector<double> brownianPath = this->Simul_Brownian(sizeSample);
//    return this->SimulPath(sizeSample,brownianPath);
//}

PoissonProcess::PoissonProcess(double intensity)
    :m_intensity(intensity)
{
}

std::vector<double> 
PoissonProcess::simulPath(const int& eventAmount) const
{
    {///simulations des temps d'un processus de Poisson simple jusqu' a un nombre d'evenements donnée
        std::vector<double> res(eventAmount, 0);
        res[0] = 0;
        for (int i = 0; i < eventAmount -1; i++)
        {
            //tau = -log(unif()) / m_intensity; ///inter-arrivé exponentielle
            res.at(i + 1) = res[i] - log(unif()) / m_intensity;// exponential
        }
        return res;
    }
}
std::vector<double>
BS_Process::simulPath(const int& sizeSample) const
{
    /// Black & Scholes log Normal paths
    /// Euler Scheme pour la discrétisation du schéma et les browniens sont obtenus avec Box-mller et incréments gaussiens
    std::vector<double> B = this->Simul_Brownian(sizeSample);
    std::vector<double> S(sizeSample);
    double dt = 1.0 / (sizeSample + 0.0);
    S[0] = getStart();
    for (int i = 1; i < sizeSample; i++)
    {
        S[i] = S[i - 1] + S[i - 1] * (getDrift() * dt + getDiffusion() * (B[i] - B[i - 1]));
    }
    return S;
}
std::vector<double> BS_Process::BS_Simul_Euler(const int& sizeSample) const
{
    /// Black & Scholes log Normal paths
    /// Euler Scheme pour la discrétisation du schéma et les browniens sont obtenus avec Box-mller et incréments gaussiens
    std::vector<double> B = this->Simul_Brownian(sizeSample);
    std::vector<double> S(sizeSample);
    double dt = 1.0 / (sizeSample + 0.0);
    S[0] = getStart();
    for (int i = 1; i < sizeSample; i++)
    {
        S[i] = S[i - 1] + S[i - 1] * (getDrift() * dt + getDiffusion() * (B[i] - B[i - 1]));
    }
    return S;
}
std::vector<double> BS_Process::BS_Simul_Euler(const int& sizeSample, const std::vector<double>& brownianPath) const
{
    /// Simulation d'une solution de l'EDS de B&S (loi log normale)
    /// Milstein scheme
    if (brownianPath.size() != sizeSample)
    {
        // If we want a simulation coarsier dont need all the point from input brownianPath
        std::vector<double> newBrownian(sizeSample);
        int inputSize = brownianPath.size();
        int startStop = 0;
        for (size_t i = 1; i < sizeSample - 1; ++i)
        {
            int startStart = startStop;
            while ((double)i / (double)sizeSample >= (double)startStop / (double)inputSize)
            {
                startStop++;
            }
            newBrownian[i] = sumSubArray(brownianPath, startStart, startStop);
        }
        newBrownian[sizeSample - 1] = sumSubArray(brownianPath, startStop, sizeSample);
        return this->fillPath(this->BS_Simul_Euler(sizeSample, newBrownian), sizeSample);
    }
    std::vector<double> S(sizeSample);
    double dt = 1.0 / (sizeSample + 0.0);
    S[0] = getStart();

    for (int i = 1; i < sizeSample; i++)
    {
        double tmp = sqrt(-2 * log(unif())) * cos(2 * PI_nb * unif()); /// N(0,1)
        S[i] = S[i - 1] + S[i - 1] * (getDrift() * dt + getDiffusion() * (brownianPath[i] - brownianPath[i - 1]));
    }
    return S;
}
std::vector<double> 
BS_Process::BS_Simul_milstein(const int& sizeSample) const
{
    std::vector<double> B = this->Simul_Brownian(sizeSample);
    return BS_Simul_milstein(sizeSample, B);
}
std::vector<double> 
BS_Process::BS_Simul_milstein(const int& sizeSample, const std::vector<double>& brownianPath) const
{
    /// Simulation d'une solution de l'EDS de B&S (loi log normale)
    /// Milstein scheme
    if (brownianPath.size() != sizeSample)
    {
        // If we want a simulation coarsier dont need all the point from input brownianPath
        std::vector<double> newBrownian(sizeSample);
        int inputSize = brownianPath.size();
        int startStop = 0;
        for (size_t i = 1; i < sizeSample - 1; ++i)
        {
            int startStart = startStop;
            while ((double)i / (double)sizeSample >= (double)startStop / (double)inputSize)
            {
                startStop++;
            }
            newBrownian[i] = sumSubArray(brownianPath, startStart, startStop);
        }
        newBrownian[sizeSample - 1] = sumSubArray(brownianPath, startStop, sizeSample);
        return this->fillPath(this->BS_Simul_milstein(sizeSample, newBrownian), sizeSample);
    }
    std::vector<double> S(sizeSample);
    double dt = 1.0 / (sizeSample + 0.0);
    S[0] = getStart();

    for (int i = 1; i < sizeSample; i++)
    {
        double tmp = sqrt(-2 * log(unif())) * cos(2 * PI_nb * unif()); /// N(0,1)
        S[i] = S[i - 1] * (1. + getDrift() * dt + getDiffusion() * (brownianPath[i] - brownianPath[i - 1])
            + (0.5 * getDiffusion() * getDiffusion() * ((brownianPath[i] - brownianPath[i - 1]) * (brownianPath[i] - brownianPath[i - 1]) - dt)));
    }
    return S;
}

void 
Simul_Sample(std::vector<double> val, int size_sample, double sd, double su, double p)
{
    int i;
    double u;
    for (i = 0; i < size_sample; i++)
    {
        u = rand() / (RAND_MAX + 1.0);
        if (u <= p)
        {
            val[i] = sd;
        }
        else
        {
            val[i] = su;
        }

    }
}

void 
display(const std::vector<double> &vect)
{
    for (size_t i = 0; i < vect.size(); i++)
    {
        std::cout << vect[i] << std::endl;
    }
}
double MAX(const double& A, const double& B)
{
    if (A < B)
        return B;
    else 
        return A;
}

double Frequence(const std::vector<double>& vect, double s)
{
    double freq = 0.;

    for (size_t i = 0; i < vect.size(); i++)
    {
        if (vect[i] == s)
        {
            freq++;
        }
    }
    freq = freq  / (vect.size() + 0.0);
    std::cout << "frequence de " << s << " = " << freq << std::endl;
    return freq;
}

std::vector<double> 
Markov_chain_3(int pathSize, int origin, double p11, double p12, double p22, double p23, double p33, double p31)
{
    // 3 States Markov chain with input transition probabilities
    std::vector<double> state(origin,pathSize);
    for (int i = 1; i < pathSize + 1; i++)
    {
        double u = rand() / (RAND_MAX + 1.0);
        if (state[i - 1] == 1)
        {//was on state 1
            if (u < p12) { state[i] = 2; }
            else if (u < p12 + p11 && u >= p12) { state[i] = 1; }
            else { state[i] = 3; }
        }
        else if (state[i - 1] == 2)
            //was on state 2
        {
            if (u < p23) { state[i] = 3; }
            else if (u < p22 + p23 && u >= p23) { state[i] = 2; }
            else { state[i] = 1; }
        }
        else
            //was on state 3
        {
            if (u < p31) { state[i] = 1; }
            else if (u < p31 + p33 && u >= p31) { state[i] = 3; }
            else { state[i] = 2; }
        }
    }
    //Frequence(state, 1);
    //Frequence(state, 2);
    //Frequence(state, 3);
    return state;
}
void Markov_chain(std::vector<std::vector<double>> p, int nb_etat, int taille, int depart)
{
    ///la voilà!!!! trajectoire markovienne en fonction d'une matrice de transition de taille quelconque
    int t;
    int i;
    int j;
    double res;
    double u;
    std::vector<double> etat(taille);
    i = 0;
    etat[0] = depart;
    std::cout << etat[0] << std::endl;
    for (t = 1; t < taille; t++) ///les instants dla chaine
    {
        i = 0;
        while (etat[t - 1] != i && i <= 10) ///parcours des etats pour trouver l'etat précédent
        {
            i++;
        }
        u = rand() / (RAND_MAX + 1.0);
        //std::cout << " u " << u << std::endl;

        res = 1; /// somme des p_i
        j = nb_etat;
        while (res >= u && j >= -2) ///dans quel intervalle de longueur P[j|i]
        {
            //std::cout << "  voila res   " << res << std::endl;
            j--;
            res = res - p[i][j];
            //std::cout <<"  res - u " <<  res - u <<" et j vo" << j << std::endl;
        }
        // std::cout <<"  voila j detat " <<  j << std::endl;
        etat[t] = j;
        //   std::cout << etat[t] << std::endl;

    }
    std::cout << "Frequences de la simulation: " << std::endl;/// derait vérifier le théoreme d'ergodicité? la loi stationnaire devrait s'afficher
    for (i = 0; i < nb_etat; i++)
    {
        Frequence(etat, i);
    }
}
std::vector<double> Simul_Brownien(int nb, double mu, double sigma) {
    std::cout << "not implemented yet" << std::endl;
    return std::vector<double>();
}
//void Gain_cumule(std::vector<double> x, int nb, double x0, double sa, double sb, double p)
//{
//    int i;
//    double u, s12;
//    x[0] = x0;
//    std::cout << x0 << std::endl;
//    for (i = 1; i < nb; i++)
//    {
//        u = rand() / (RAND_MAX + 1.0);
//        s12 = (u <= 1 - p) ? sa : sb;
//        x[i] = x[i - 1] + s12;
//        std::cout << x[i] << std::endl;
//    }
//}
//void M_C_Gain_cumule(std::vector<double> x, int nb, int nb_mc, double x0, double sa, double sb, double p)
//{
//    /// Exo de Deauville.
//    int i;
//    int j;
//    double u, s12;
//    std::vector<double> t;
//    t = new double[nb_mc];
//    x[0] = x0;
//    std::cout << x0 << std::endl;
//    for (j = 0; j < nb_mc; j++)
//    {
//        for (i = 1; i < nb; i++)
//        {
//            u = rand() / (RAND_MAX + 1.0);
//            s12 = (u <= 1 - p) ? sa : sb;
//            x[i] = x[i - 1] + s12;
//            //std::cout << x[i] << std::endl ;
//        }
//        i = 1;
//        t[j] = x[nb - 1];
//        std::cout << t[j] << std::endl;
//    }
//    delete[]t;
//}
//void M_C_Gain_cumule_Max(std::vector<double> x, int nb, int nb_mc, int maxx, double x0, double sa, double sb, double p)
//{///sensé afficher nb la proba du max
//    int i;
//    int j;
//    double u, s12;
//    std::vector<double> t;
//    t = new double[nb_mc];
//    x[0] = x0;
//    std::cout << x0 << std::endl;
//    for (j = 0; j < nb_mc; j++)
//    {
//        for (i = 1; i < nb; i++)
//        {
//            u = rand() / (RAND_MAX + 1.0);
//            s12 = (u <= 1 - p) ? sa : sb;
//            if (x[i - 1] <= maxx)
//            {
//                x[i] = x[i - 1] + s12;
//            }
//            else
//            {
//                x[i] = maxx;//std::cout << x[i] << std::endl ;
//            }
//            i = 1;
//            t[j] = x[nb - 1];
//            std::cout << t[j] << std::endl;
//        }
//        delete[]t;
//    }
//}
//
//void Fcts_Barriere_up_and_down(std::vector<double> val, int taille, double sup, double inf)
//{ ///fonction suspecte
//    int i;
//    for (i = 0; i < taille; i++)
//    {
//        if (val[i] <= inf) { val[i] = inf; }
//        else if (val[i] >= sup) { val[i] = sup; }
//        std::cout << val[i] << std::endl;
//
//    }
//}
//
//std::vector<double> Temps_Atteinte(int nb, int nb_stop, double x0, double sa, double sb, double p)
//{
//    /// simulation d'un temps d'atteinte d'une marche aléatoire assymétrique
//    int k, j;
//    j = 0;
//    double u, res, tau, s12;
//    tau = 0;
//    std::vector<double> t;
//    t = new double[nb];
//    for (k = 0; k < nb; k++)
//    {
//        while (tau != 3 && j < nb_stop)
//        {
//            u = rand() / (RAND_MAX + 1.0);
//            s12 = (u <= 1 - p) ? sa : sb;
//            tau += s12;
//            j++;
//            std::cout << "tau 11 " << tau << std::endl;
//
//            //        std::cout << j << std::endl;
//            //        std::cout << tau << std::endl;
//        }
//        //    std::cout << "j " << j << std::endl;
//        if (j < nb_stop)
//        {
//            t[k] = j;
//        }
//        else
//        {
//            t[k] = -1;
//        }
//        std::cout << tau << std::endl;
//        tau = 0;
//        j = 0;
//    }
//    return t;
//}
std::vector<double> Proc_O_U_Simul(int nb, double origine, double moyenne, double diffusion, double speed_retour)
{
    /// Simulation d'un processus Orstein-Uhlenbeck
    /// Bon pour observer le retour à la moyenne, faire du monte-carlo sur les taux
    std::vector<double> B(nb);
    std::vector<double> X(nb);
    B = Simul_Brownien(nb, 0, 1);
    X[0] = origine;
    std::cout << X[0] << std::endl;
    for (int i = 1; i < nb; i++)
    {
        X[i] = X[i - 1] - speed_retour * (X[i - 1] - moyenne) * (1.0 / (nb + 0.0)) + diffusion * (B[i] - B[i - 1]);
        std::cout << X[i] << std::endl;
    }
    return X;
}
std::vector<double> Proc_CIR_Simul(int nb, double origine, double moyenne, double diffusion, double speed_retour)
{   /// Simulation d'un processus CIR
    /// Bon pour observer le retour la positivité des taux ssi sigma^2<2*k*theta, n'hésite pas a faire du monte-carlo en moyennant les sorties de ce truc
        /// on récupère en sortie une trajectoire de CIR

    std::vector<double> B(nb);
    std::vector<double> X(nb);
    B = Simul_Brownien(nb, 0, 1);
    X[0] = origine;
    std::cout << X[0] << std::endl;
    for (int i = 1; i < nb; i++)
    {
        X[i] = X[i - 1] - speed_retour * (X[i - 1] - moyenne) * (1.0 / (nb + 0.0)) + diffusion * sqrt(X[i]) * (B[i] - B[i - 1]);
        std::cout << X[i] << std::endl;
    }
    return X;
}
std::vector<double> Proc_pseudo_CEV(int nb, double origine, double moyenne, double nu, double puissance)
{
    ///Petit modèle de vol locale qui est équivalent à un CEV pour de grandes valeurs de sous jacent mais évite le probleme de la division par zéro quand  S est petit

    /// Par une discrétisation Eulérienne; les brownien sont obtenus par incréments gaussiens et Box Muller

    /// on récupère en sortie une trajectoire de pseudo-CEV
    std::vector<double> B(nb);
    std::vector<double> X(nb);
    B = Simul_Brownien(nb, 0, 1);
    X[0] = origine;
    std::cout << X[0] << std::endl;
    for (int i = 1; i < nb; i++)
    {//     res = sqrt(dt)*Nu*(x.*x.^(Delta))./(sqrt(1+x.^2));
        X[i] = X[i - 1] + moyenne * (X[i - 1]) * (1.0 / (nb + 0.0)) + (B[i] - B[i - 1]) * (sqrt((1.0 / (nb + 0.0)))) * nu * pow(X[i - 1], puissance + 1) / (sqrt(1 + X[i - 1]));
        std::cout << X[i] << std::endl;
    }
    return X;
}
double sd_hat(std::vector<double> X, int n)
{
    /// Calcul de l'estimateur empirique de la variance

    double tmp_moyenne(0);
    double res(0);
    //  std::cout << " debut de std " << res << std::endl;

    for (int i = 0; i < n; i++)
    {
        tmp_moyenne += X[i];
    }
    tmp_moyenne = tmp_moyenne / (double)n;
    for (int i = 0; i < n; i++)
    {
        res = res + pow(X[i] - tmp_moyenne, 2);
    }
    res = res / (n - 1.0);
    res = sqrt(res);
    return res;
}
std::vector<double> Trier(std::vector<double> tab, int taille)
{
    /// fonction de tri non optimale but still efficace
    int i;
    for (i = 0; i < taille; i++)
    {
        int j;
        double x;
        x = tab[i];
        j = i;

        while ((j > 0) && tab[j - 1] > x)
        {
            tab[j] = tab[j - 1];
            j--;
        }
        tab[j] = x;
    }
    return tab;
}
int TVI(std::vector<double> t, int taille, double e)
{
    /// sort l'indice ou on a la valeur la plus proche de l'argument en entrée
    /// on va trier le tableau anyway donc n'importe quel tableau en entrée fera lafaire
    int inf(0);
    int sup(taille - 1);
    std::vector<double> sort_tab;
    sort_tab = Trier(t, taille);
    int tmp;
    tmp = taille / 2;
    while (/*inf<sup ||*/ abs(sort_tab[tmp] - e) > 0.5)
    {
        //std::cout << " dans le while " <<std::endl;

        if (e <= sort_tab[tmp])
        {
            sup = tmp;
            tmp = (tmp - inf + 1) / 2;
        }
        else
        {
            inf = tmp;
            tmp = (sup - tmp + 1) / 2;
        }

    }
    // std::cout << "la plus proche valeur de " << e << " est " << sort_tab[tmp] << " atteint en " << tmp << std::endl;
    return tmp;
}
//double * Call_PCEV_MC(int nb_MC,int nb, double origine, double moyenne, double nu, double puissance)
//{
//double *MC = new double [nb_MC];
//   double *a;
//    for(int i(0);i<nb_MC;i++)
//    {
//       a = Proc_pseudo_CEV(int nb, double origine, double moyenne, double nu, double puissance);
//       MC[i] = -1;
//    }
//
//}


double f(double x, double lambda)
{
    double res;
    res = exp(x * lambda);
    return res;
}
double h(double x)
{
    double res;
    res = 5 * sqrt(1 + x);
    return res;
}

//void Call_BS_MC(int nb_MC, int nb_pas, double spot, double mu, double diffusion, double K)
//{
//    ///affiche le prix par M-C d'un call sous BS par schéma d'Euler
//    std::vector<double> MC(nb_MC);
//    std::vector<double> S(nb_pas);
//    double tmp;
//    double tmpp;
//    double tmpm;
//    tmp = 0;
//    int i;
//    for (i = 0; i < nb_MC; i++)
//    {
//        S = BS_Simul(nb_pas, spot, mu, diffusion);
//        MC[i] = S[nb_pas - 1];
//        std::cout << " MC[i] " << S[nb_pas - 1] << std::endl;
//        //tmp = tmp + S[nb_pas - 1]/nb_MC;
//    }         // std::cout << " valeur initiale payoff  " << MAX(MC[0] - K,0) << std::endl;
//
//    for (i = 0; i < nb_MC; i++)///On fait la moyenne des S qui doit converger
//    {
//        tmp = tmp + MAX(MC[i] - K, 0) / (nb_MC + 1.0 - 1.0);/// Pay-off d'un call = (St-K)^+
//        std::cout << " tmp " << tmp << std::endl;
//    }
//    //tmp = tmp/nb_MC;
//    tmpp = tmp + 1.96 * sd_hat(MC, nb_MC) / sqrt(nb_MC);
//    tmpm = tmp - 1.96 * sd_hat(MC, nb_MC) / sqrt(nb_MC);
//
//    tmp = exp(-mu) * tmp;
//    std::cout << " prix bs " << tmp << std::endl;
//    std::cout << " Intervalle de confiance a 5% : " << tmpm << " ; " << tmpp << std::endl;
//}


