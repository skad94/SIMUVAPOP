#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>
#include <vector>
#include "Simulations_VA.h"
#define PI_nb 3.14159

double * Simulation (int size_sample)
{
    /// génére un tableau de taille donnée
    double *t;
    t = new double[size_sample];
    return t;
}

void Free_Simulation(double *t)
{
    /// à utiliser après une simulation pour libérer
    delete[] t;
}

void Initialise(double *vect, int size_sample, double val=1.1)
{
    int i;
    for(i=0; i<size_sample ;i++)
    {
        vect[i]=val;
    }
}


void Simul_Sample (double * val, int size_sample, double sd, double su, double p)
{
    int i;
    double u;
    for(i=0;i<size_sample;i++)
    {
        u=rand()/(RAND_MAX +1.0);
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

void affiche(double * vect, int size_sample)
{
    /// exo de mon directeur
    int i;
    for(i=0;i<size_sample;i++)
    {
        std::cout << vect[i] << std::endl;
    }
}
double MAX(double A, double B)
{
 double res;
 if (A<B){res = B;}
 else {res = A;}
 return res;
}

double Frequence(double *vect, int size_sample, double s)
{
    int i;
    double freq = 0;

    for(i=0 ; i < size_sample ; i++)
    {
        if (vect[i]==s)
        {
          freq ++;
        }
       // freq = freq + ((vect[i]==s)? 1.0f : 0.0f)/ size_sample;
    }
freq = (freq+0.0)/(size_sample +0.0);
std::cout << "frequence de " << s << " = "<<freq << std::endl;
    return freq;
}

double * Proc_Poisson_Simul(double intensite, int nb_event)
 {///simulations des temps d'un processus de Poisson simple jusqu' a un nombre d'evenements donnée
  /// ATTENTION  a la libération de l'espace memoire
    int i;
    double tau;
    double *t;
    t = new double[nb_event +1];
    t[0]=0;
    for(i=0;i<nb_event +1;i++)
    {
       tau = rand()/(RAND_MAX +1.0);
       tau = -1*log(tau)/intensite; ///inter-arrivé exponentielle
       t[i+1] = t[i] + tau;
       std::cout << t[i] << std::endl;
    }
    return t;
 }

double * Simul_Brownien(int nb,double mu, double sigma)
{
    /// Simulation d'un brownien par les incréments gaussiens indépendants
    int i;
    double *n;
    n = new double[nb];
    double *W;
    W = new double[nb];
    double x1;
    double x2;
    x1=0;
    x2=0;

    /// Simulation d'un échantillon 1D normale centré réduit par box Muller
    for(i=0;i<nb;i++)
    {
        x1 = rand()/(RAND_MAX +1.0);
        x2 = rand()/(RAND_MAX +1.0);
        n[i] = sqrt(-2*log(x1))*cos(2*PI_nb*x2);
        n[i] = sigma*n[i] + mu;
        //std::cout << n[i] << std::endl;
    }

    W[0]=0;
   // std::cout << W[0] << std::endl;
    for(i=1;i<nb;i++)
    {
        W[i] = W[i-1] + sqrt((1.0+0.0)/(nb+1.1-1.1))*n[i-1];
       // std::cout << W[i] << std::endl;
    }
    delete[] n;
    return W;
}

void Markov_chain_3(int nb, int depart,double p11, double p12, double p22, double p23, double p33, double p31)
{
    ///Chaine de markov à 3 états (en attendant la version généralisé...to be continued) qui prend les proba de transition et affiche une trajectoire
    int i;
    double u;
    double *etat;
    etat = new double [nb];
    etat[0] = depart;
    for(i=1;i<nb+1;i++)
    {
        u=rand()/(RAND_MAX +1.0);
       std::cout << i << ") ";

        if(etat[i-1]==1)
        {///on etait a letat 1
            if(u<p12){etat[i]=2;}
            else if(u<p12+p11 && u>=p12){etat[i]=1;}
            else {etat[i]=3;}
        }
        else if (etat[i-1]==2)
        ///on etait a letat 2
        {
            if(u<p23){etat[i]=3;}
            else if(u<p22+p23 && u>=p23){etat[i]=2;}
            else {etat[i]=1;}
        }
        else
        ///on etait a letat 3
        {
            if(u<p31){etat[i]=1;}
            else if(u<p31+p33 && u>=p31){etat[i]=3;}
            else {etat[i]=2;}
        }
       std::cout << etat[i] << std::endl;
    }
    Frequence(etat,nb,1);
    Frequence(etat,nb,2);
    Frequence(etat,nb,3);
     delete[] etat;
}
void Markov_chain(double **p, int nb_etat, int taille, int depart)
{
    ///la voilà!!!! trajectoire markovienne en fonction d'une matrice de transition de taille quelconque
    int t;
    int i;
    int j;
    double res;
    double u;
    double *etat;
    i=0;
    etat = new double [taille];
    etat[0] = depart;
    std::cout << etat[0] << std::endl;
    for(t=1;t<taille;t++) ///les instants dla chaine
    {
        //std::cout << i <<"  voila i debut de for  " << std::endl;
        i=0;
        //std::cout << i <<"  voila i ki vo zero  " << std::endl;
       while(etat[t-1]!=i && i <= 10) ///parcours des etats pour trouver l'etat précédent
        {
            i++;
            //    std::cout << i <<"  voila i du while  " << std::endl;

        }
//        if(i==10)
//        {
//                std::cout << "probleme" << std::endl;
//                delete[] etat;
//
//        }
//        else
//        {
            u=rand()/(RAND_MAX+1.0);
            //std::cout << " u " << u << std::endl;

        res=1; /// somme des p_i
        j=nb_etat;
       while(res>=u && j>=-2) ///dans quel intervalle de longueur P[j|i]
       {
           //std::cout << "  voila res   " << res << std::endl;
           j--;
           res = res - p[i][j];
           //std::cout <<"  res - u " <<  res - u <<" et j vo" << j << std::endl;
       }
          // std::cout <<"  voila j detat " <<  j << std::endl;
       etat[t]=j;
    //   std::cout << etat[t] << std::endl;

        }
   std::cout << "Frequences de la simulation: " << std::endl;/// derait vérifier le théoreme d'ergodicité? la loi stationnaire devrait s'afficher
  for(i=0;i<nb_etat;i++)
  {
      Frequence(etat,taille,i);
  }
delete[] etat;
    }
    //}

void Gain_cumule(double *x, int nb, double x0, double sa, double sb, double p)
{
    int i;
    double u, s12;
    x[0]=x0;
    std::cout << x0 << std::endl;
    for(i=1;i<nb;i++)
    {
        u=rand()/(RAND_MAX+1.0);
        s12 = (u<=1-p)?sa:sb;
        x[i]=x[i-1]+s12;
        std::cout << x[i] << std::endl;
    }
}
void M_C_Gain_cumule(double *x, int nb,int nb_mc, double x0, double sa, double sb, double p)
{
    /// Exo de Deauville.
    int i;
    int j;
    double u, s12;
    double *t;
    t = new double[nb_mc];
    x[0]=x0;
    std::cout << x0 << std::endl;
    for(j=0;j<nb_mc;j++)
    {
        for(i=1;i<nb;i++)
    {
        u=rand()/(RAND_MAX+1.0);
        s12 = (u<=1-p)?sa:sb;
        x[i]=x[i-1]+s12;
        //std::cout << x[i] << std::endl ;
    }
    i=1;
    t[j] = x[nb-1];
    std::cout << t[j] << std::endl;
    }
    delete []t;
}
void M_C_Gain_cumule_Max(double *x, int nb,int nb_mc,int maxx, double x0, double sa, double sb, double p)
{///sensé afficher nb la proba du max
    int i;
    int j;
    double u, s12;
    double *t;
    t = new double[nb_mc];
    x[0]=x0;
    std::cout << x0 << std::endl;
    for(j=0;j<nb_mc;j++)
    {
        for(i=1;i<nb;i++)
    {
        u=rand()/(RAND_MAX+1.0);
        s12 = (u<=1-p)?sa:sb;
        if(x[i-1] <= maxx)
        {
        x[i]=x[i-1]+s12;
        }
        else
          {
              x[i] = maxx;//std::cout << x[i] << std::endl ;
          }
    i=1;
    t[j] = x[nb-1];
    std::cout << t[j] << std::endl;
    }
    delete []t;
}
}

void Fcts_Barriere_up_and_down(double * val, int taille, double sup, double inf)
{ ///fonction suspecte
    int i;
for(i=0;i<taille;i++)
   {
      if(val[i] <= inf) {val[i]=inf;}
      else if(val[i]>= sup) {val[i]=sup;}
    std::cout << val[i] << std::endl;

   }
}

double * Temps_Atteinte(int nb,int nb_stop, double x0, double sa, double sb, double p)
{
    /// simulation d'un temps d'atteinte d'une marche aléatoire assymétrique
    int k,j;
    j=0;
    double u, res, tau, s12;
    tau = 0;
    double *t;
    t = new double [nb];
    for(k=0;k<nb;k++)
    {
    while(tau!=3 && j<nb_stop)
    {
        u=rand()/(RAND_MAX+1.0);
        s12 = (u<=1-p)?sa:sb;
        tau += s12;
        j++;
        std::cout << "tau 11 " << tau << std::endl;

//        std::cout << j << std::endl;
//        std::cout << tau << std::endl;
    }
//    std::cout << "j " << j << std::endl;
    if(j<nb_stop)
    {t[k] = j;}
    else
        {t[k] = -1;}
    std::cout << tau << std::endl;
    tau=0;
    j=0;
    }
    return t;
}
double * Proc_O_U_Simul(int nb, double origine, double moyenne, double diffusion, double speed_retour)
{
    /// Simulation d'un processus Orstein-Uhlenbeck
    /// Bon pour observer le retour à la moyenne, faire du monte-carlo sur les taux
    double * B;
    double * X;
    int i;
    X = new double [nb];
    B = Simul_Brownien(nb,0,1);
    X[0]= origine;
    std::cout << X[0] << std::endl;
    for(i=1;i<nb;i++)
      {
       X[i] = X[i-1] - speed_retour * (X[i-1]-moyenne)*(1.0/(nb+0.0)) + diffusion* (B[i]-B[i-1]);
       std::cout << X[i] << std::endl;
      }
    Free_Simulation(B);
    return X;
}
double * Proc_CIR_Simul(int nb, double origine, double moyenne, double diffusion, double speed_retour)
{   /// Simulation d'un processus CIR
    /// Bon pour observer le retour la positivité des taux ssi sigma^2<2*k*theta, n'hésite pas a faire du monte-carlo en moyennant les sorties de ce truc
        /// on récupère en sortie une trajectoire de CIR

    double * B;
    double * X;
    int i;
    X = new double [nb];
    B = Simul_Brownien(nb,0,1);
    X[0]= origine;
    std::cout << X[0] << std::endl;
    for(i=1;i<nb;i++)
      {
       X[i] = X[i-1] - speed_retour * (X[i-1]-moyenne)*(1.0/(nb+0.0)) + diffusion* sqrt(X[i])*(B[i]-B[i-1]);
       std::cout << X[i] << std::endl;
      }
    Free_Simulation(B);
    return X;
}
double * Proc_pseudo_CEV(int nb, double origine, double moyenne, double nu, double puissance)
{
    ///Petit modèle de vol locale qui est équivalent à un CEV pour de grandes valeurs de sous jacent mais évite le probleme de la division par zéro quand  S est petit

    /// Par une discrétisation Eulérienne; les brownien sont obtenus par incréments gaussiens et Box Muller

    /// on récupère en sortie une trajectoire de pseudo-CEV
    double * B;
    double * X;
    int i;
    X = new double [nb];
    B = Simul_Brownien(nb,0,1);/// Brownien par incréments gaussiens et Box Muller
    X[0]= origine;
    std::cout << X[0] << std::endl;
    for(i=1;i<nb;i++)
      {//     res = sqrt(dt)*Nu*(x.*x.^(Delta))./(sqrt(1+x.^2));
       X[i] = X[i-1] + moyenne*(X[i-1])*(1.0/(nb+0.0)) + (B[i]-B[i-1])*(sqrt((1.0/(nb+0.0))))*nu*powf(X[i-1],puissance+1)/(sqrt(1 + X[i-1]));
       std::cout << X[i] << std::endl;
      }
    Free_Simulation(B);
    return X;
}
double sd_hat(double *X,int n)
{
    /// Calcul de l'estimateur empirique de la variance

    double tmp_moyenne(0);
    double tmp_carre;
    double res(0);
          //  std::cout << " debut de std " << res << std::endl;

    for(int i=0;i<n;i++)
    {
        tmp_moyenne +=  X[i];
    }
      tmp_moyenne = tmp_moyenne/n;
       for(int i=0;i<n;i++)
      {
        res = res + pow(X[i] - tmp_moyenne,2);
      }
    res = res/(n-1);
    res = sqrt(res);
    return res;
}
double * Trier(double *tab, int taille)
{
    /// fonction de tri non optimale but still efficace
    int i;
    for (i=0;i<taille;i++)
        {
        int j;
        double x;
        x = tab[i];
        j = i;

        while((j>0) && tab[j-1]>x)
            {
            tab[j] = tab[j-1];
            j--;
            }
        tab [j] = x;
        }
    return tab;
}
int TVI(double *t,int taille,double e)
{
    /// sort l'indice ou on a la valeur la plus proche de l'argument en entrée
    /// on va trier le tableau anyway donc n'importe quel tableau en entrée fera lafaire
    int inf(0);
    int sup(taille - 1);
    double * sort_tab;
    sort_tab = Trier(t,taille);
    int tmp;
    tmp = taille/2;
    while (/*inf<sup ||*/ abs(sort_tab[tmp]-e)>0.5)
    {
        //std::cout << " dans le while " <<std::endl;

        if(e<=sort_tab[tmp])
            {
                sup = tmp;
                tmp = (tmp - inf +1)/2;
            }
        else
            {
                inf = tmp;
                tmp = (sup - tmp +1)/2;
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

double * BS_Simul(int nb, double origine, double mu, double diffusion)
{
    /// Simulation d'une trajectoire de l'EDS de B&S (loi log normale)
    /// Schéma d'Euler pour la discrétisation du schéma et les browniens sont obtenus avec box-mller et incréments gaussiens
    double * B;
    double * S;
    int i;
    double dt;
    dt = 1.0/(nb+0.0);
    S = new double [nb];
    B = Simul_Brownien(nb,0,1);
    S[0]= origine;
   // std::cout << S[0] << std::endl;
    for(i=1;i<nb;i++)
      {
       S[i] = S[i-1] + S[i-1]* (mu*dt + diffusion*(B[i]-B[i-1]));
       //std::cout << S[i] << std::endl;
      }
    Free_Simulation(B);
    return S;
}
double * BS_Simul_milstein(int nb, double origine, double mu, double diffusion)///PAS FINIE DU TOUT
{
    /// Simulation d'une solution de l'EDS de B&S (loi log normale)
    /// Schéma de Milstein pour la discrétisation du schéma et les browniens sont obtenus avec box-muller et incréments gaussiens
    double * B;
    double * S;
    int i;
    double x1;
    double x2;
    double dt;
    double tmp(0);
    dt = 1.0/(nb+0.0);
    S = new double [nb];
    B = Simul_Brownien(nb,0,1);
    S[0]= origine;

    /// Simulation d'un échantillon 1D normale centré réduit par box Muller

    for(i=1;i<nb;i++)
      {
        x1 = rand()/(RAND_MAX +1.0);
        x2 = rand()/(RAND_MAX +1.0);
        tmp = sqrt(-2*log(x1))*cos(2*PI_nb*x2); /// tirage N(0,1)
        S[i] = S[i-1];//( 1 + diffusion*sqrt(dt)*tmp + ( pow(tmp,2) - 1)*(pow(diffusion,2)*dt)/(2.0));/// ajout de l'innovation du crochet
       //std::cout << S[i] << std::endl;
      }
    Free_Simulation(B);
    return S;
}
double f(double x,double lambda)
{
    double res;
    res = exp(x*lambda);
    return res;
}
double h(double x)
{
    double res;
    res = 5*sqrt(1+x);
    return res;
}
double * Lyvat_Simul(int nb, double origine, double mu, double diffusion)
{
    /// Simulation des differents trucs de lyvath
    /// double * Proc_Poisson_Simul(double intensite, int nb_event)
///double * BS_Simul(int nb, double origine, double mu, double diffusion)

    double * T;
//    double * S;
//    double * Q;
//    double * B;
//    double * Z;
//    int i;
//    double dt;
//    dt = 1.0/(nb+0.0);
//    S = new double [nb];
//    T = Proc_Poisson_Simul(intensite,nb_event);
//    B = Simul_Brownien(nb,0,1);
//    S = BS_Simul(int nb, double origine, double mu, double diffusion);
//
//    /Q ou le quantity de dividends
//    for(i=0;i<nb;i++)
//    {
//        Q[i] = q;
//    }
//    /Z ou le quantity de dividends over le tps
//    Z[0] = q;
//    for(i=0;i<nb;i++)
//    {
//        Z[i+1] = Z[i]+q;
//    }
//    /X ou le cash de l'entreprise
//    X[0] = 0;
//      for(i=0;i<nb;i++)
//    {
//        X[i+1] = r*X[i]*dt + h(Q[i])(b*dt + eta *B[i+1]-B[i]) - (Z[i+1] - Z[i]) ;
//    }
//    Free_Simulation(B);
    return T;
}
double * Lyvat_Simul_Q(int nb, double origine, double mu, double diffusion)
{
    /// Simulation d'une trajectoire de l'EDS de B&S (loi log normale)
    /// Schéma d'Euler pour la discrétisation du schéma et les browniens sont obtenus avec box-mller et incréments gaussiens
    double * B;
    double * S;
    int i;
    double dt;
    dt = 1.0/(nb+0.0);
    S = new double [nb];
    B = Simul_Brownien(nb,0,1);
    S[0]= origine;
   // std::cout << S[0] << std::endl;
    for(i=1;i<nb;i++)
      {
       S[i] = S[i-1] + S[i-1]* (mu*dt + diffusion*(B[i]-B[i-1]));
       //std::cout << S[i] << std::endl;
      }
    Free_Simulation(B);
    return S;
}

double * Lyvat_Simul_X(int nb, double origine, double mu, double diffusion)
{
    /// Simulation d'une trajectoire de l'EDS de B&S (loi log normale)
    /// Schéma d'Euler pour la discrétisation du schéma et les browniens sont obtenus avec box-mller et incréments gaussiens
    double * B;
    double * S;
    int i;
    double dt;
    dt = 1.0/(nb+0.0);
    S = new double [nb];
    B = Simul_Brownien(nb,0,1);
    S[0]= origine;
   // std::cout << S[0] << std::endl;
    for(i=1;i<nb;i++)
      {
       S[i] = S[i-1] + S[i-1]* (mu*dt + diffusion*(B[i]-B[i-1]));
       //std::cout << S[i] << std::endl;
      }
    Free_Simulation(B);
    return S;
}
void Call_BS_MC(int nb_MC, int nb_pas, double spot, double mu, double diffusion, double K)
{
    ///affiche le prix par M-C d'un call sous BS par schéma d'Euler
     double *MC;
        MC = new double [nb_MC];
        double *S;
        S = new double [nb_pas];
        double tmp ;
        double tmpp ;
        double tmpm ;
        tmp = 0;
        int i;
    for(i=0;i<nb_MC;i++)
    {
     S = BS_Simul(nb_pas,spot,mu,diffusion);
     MC[i] = S[nb_pas - 1];
     std::cout << " MC[i] " << S[nb_pas - 1] << std::endl;
     //tmp = tmp + S[nb_pas - 1]/nb_MC;
     }         // std::cout << " valeur initiale payoff  " << MAX(MC[0] - K,0) << std::endl;

      for(i=0;i<nb_MC;i++)///On fait la moyenne des S qui doit converger
      {
          tmp = tmp + MAX(MC[i] - K,0)/(nb_MC + 1.0 - 1.0);/// Pay-off d'un call = (St-K)^+
          std::cout << " tmp " << tmp << std::endl;
      }
        //tmp = tmp/nb_MC;
        tmpp = tmp + 1.96*sd_hat(MC,nb_MC)/sqrt(nb_MC);
        tmpm = tmp - 1.96*sd_hat(MC,nb_MC)/sqrt(nb_MC);

        Free_Simulation(S);
           tmp = exp(-mu)*tmp;
           std::cout << " prix bs " << tmp << std::endl;
           std::cout << " Intervalle de confiance a 5% : " << tmpm << " ; " << tmpp << std::endl;

                   Free_Simulation(MC);

}
