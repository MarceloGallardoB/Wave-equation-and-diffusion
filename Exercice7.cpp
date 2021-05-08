// ============================================================================================================================================================
// Import des bibliothèques et création de type
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.h"

using namespace std;
typedef enum{not_defined,dirichlet,neumann,harmonic} bound_cond;
#define M_PI       3.14159265358979323846


// ============================================================================================================================================================
// Méthodes générales
// ----------------------------------------------------------------------------------------------------------------------
// Méthode pour lire les conditions aux limites
bound_cond read_boundary_condition(string const& bnd_in){
  if(bnd_in == "dirichlet")
    return dirichlet;
  else if(bnd_in == "neumann")
    return neumann;
  else if(bnd_in == "harmonic")
    return harmonic;
  else{
    cerr << "Please select "+bnd_in+"=""dirichlet"", ""neumann"", ""harmonic""." << endl;
    return not_defined;
  }
}

// ----------------------------------------------------------------------------------------------------------------------
// TODO : Calcul de l'énergie de l'onde
double energy(vector<vector<double>> const& f, double const& dx, double const& dy,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD){
  double energy_(0.);
  // TODO: Compléter ici en utilisant la formule des trapèzes
  for(double x(start_xiD); x < end_xiD ; x += 1){
    for(double y(start_yiD); y < end_yiD ; y += 1){
      energy_ += (dx*dy/4)*(pow(f[x][y], 2) + pow(f[x+1][y], 2) + pow(f[x][y+1], 2) + pow(f[x+1][y+1], 2));
    }
  }

  return energy_;
}

// ----------------------------------------------------------------------------------------------------------------------
// TODO: Calculer la perturbation
double perturbation(double const& t,double const& x,double const& y,\
double const& xL, double const& xR,double const& yL, double const& yU, \
double const& pert_velocity, double const& pert_amplitude, \
int const& mode_num_x,int const& mode_num_y, bool const& f_tilde){
  if (f_tilde){
    return pert_amplitude*(cos((mode_num_x*M_PI*x/(xR - xL)) - (mode_num_y*M_PI*y/(yU - yL))) - cos((mode_num_x*M_PI*x/(xR - xL)) + (mode_num_y*M_PI*y/(yU - yL))));
  }
  else{
    return pert_amplitude*sin(pert_velocity*t);
  }
}

// ----------------------------------------------------------------------------------------------------------------------
// Surcharge de l'opérateur pour écrire les élements d'un tableau 1D
template <class T> ostream& operator<< (ostream& o, vector<T> const& v){
  unsigned int len(v.size());

  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";

  if(len > 0)
    o << v[len-1];

  return o;
}

// ----------------------------------------------------------------------------------------------------------------------
// Surcharge de l'opérateur pour écrire les élements d'un tableau 2D
template<typename T> void print_table(ostream& o, double const& t , vector<vector<T>> const& v,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD){
  unsigned int i,j;
  
  for(i = start_xiD; i < end_xiD + 1; ++i){
      o << t << " ";
      for(j = start_yiD; j < end_yiD; ++j){
        o << v[i][j] << " "; 
      }
      o << v[i][end_yiD] << endl;
  }
}



// ============================================================================================================================================================
// Fonction u^2(x,y)
// ----------------------------------------------------------------------------------------------------------------------
// Classe abstraite
class U2 {
public:
  U2(ConfigFile const& configFile){
    xL=configFile.get<double>("xL"); 
    xR=configFile.get<double>("xR");
    yL=configFile.get<double>("yL");
    yU=configFile.get<double>("yU");
  }

  double get_left_extremum(){
    return xL;
  }
  double get_right_extremum(){
    return xR;
  }
  double get_upper_extremum(){
    return yU;
  }
  double get_lower_extremum(){
    return yL;
  }
  // Méthodes virtuelles pures => classe abstraite
  virtual double operator()(double const& x, double const& y) const = 0; // Evalue u^2 au point (x,y)
  virtual double max() const = 0; // Renvoie max(u^2(x,y))

protected:
  double xL,xR,yL,yU;
};

// ----------------------------------------------------------------------------------------------------------------------
// Fonction constante
class U2_const: public U2 {
public:
  // Pas de constructeur par défaut => on force à spécifier une valeur
  U2_const(ConfigFile const& configFile) :
  U2(configFile), u2(pow(configFile.get<double>("u"),2)){}

  // Définition des méthodes virtuelles pures :
  double operator()(double const& x,double const& y) const{
    return u2;
  }

  double max() const{
    return u2;
  }

private:
  double u2;
};

// ----------------------------------------------------------------------------------------------------------------------
// Fonction sous forme d'onde
class U2_onde: public U2 {
public:
  U2_onde(ConfigFile const& configFile):
  U2(configFile),
  g(configFile.get<double>("g")),
  h0(configFile.get<double>("h0")),
  h1(configFile.get<double>("h1")),
  a(configFile.get<double>("a")),
  b(configFile.get<double>("b")){}

protected:
  double g, h0, h1, a, b, Lx;
};

// ----------------------------------------------------------------------------------------------------------------------
// Onde cas 1
class U2_onde_cas1: public U2_onde {
public:
  U2_onde_cas1(ConfigFile const& configFile):
  U2_onde(configFile) {}
  
  double h(double const& x,double const& y) const{
    double h_val(h0);
    if(x>a && x<b) h_val = h0+(h1-h0)*sin(M_PI*(x-a)/(b-a));
      return h_val;
  }
  
  double h_prime(double const& x) const{
    return (M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a)))/(b-a);
  }
  
  double h_prime2(double const& x) const{
    return -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a)))/((b-a)*(b-a));
  }
  
  double operator()(double const& x,double const& y) const{
    return g * h(x,y);
  }
  
  double max() const{
    unsigned int iteration(0),maxit(1000);
    double error(1.e20),tol(1.e-10);
    double x(5.e-1*(a+b));

   // Use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){
     // Compute new solution
     x = x - (h_prime(x)/h_prime2(x));
     // Compute new error
     error = abs(h_prime(x));
     // Increase iterator
     iteration += 1;
   }

   // Send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // Check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;
   } 
   else {
     if(h_prime2(x)<0.e0){
       cout << "Newton method: maximum not found!" << endl;
     }
   }
   // Calculer la vitesse maximale
   return g*(std::max(h0,std::max(h(std::max(a,std::min(x,b)),0.),h1)));
  }
};

// ----------------------------------------------------------------------------------------------------------------------
// Onde cas 2
class U2_onde_cas2: public U2_onde {
protected:
  double Ly;
public:
  U2_onde_cas2(ConfigFile const& configFile) :
  U2_onde(configFile), Ly(configFile.get<double>("Ly")) {}
  
  double h(double const& x,double const& y) const {
    double h_value(h0);
    if(x>a && x<b ) h_value = h0+(h1-h0)*sin(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly);
      return h_value;
  }

  vector<double> h_prime(double const& x,double const& y) const{
    vector<double> grad(2);
    grad[0] = (M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/(b-a);
    grad[1] = (M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*cos(M_PI*y/Ly))/Ly;
    return grad;
  }

  vector<double> h_prime2(double const& x,double const& y) const{
    vector<double> d2h(3,0.0);
    d2h[0] = -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/((b-a)*(b-a));
    d2h[1] = (M_PI*M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a))*cos(M_PI*y/Ly))/((b-a)*Ly);
    d2h[2] = -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/(Ly*Ly);
    return d2h;
  }

  vector<vector<double>> h_prime2_inv(double const& x,double const& y) const{
    vector<vector<double>> jac_inv(2,vector<double>(2));
    vector<double> d2h(h_prime2(x,y));
    double den=0.0;
    
    // Calculer l'inverse du jacobien
    den = d2h[0]*d2h[2]-d2h[1]*d2h[1];
    jac_inv[0][0] = d2h[2]/den;
    jac_inv[0][1] = -d2h[1]/den;
    jac_inv[1][0] = d2h[1];
    jac_inv[1][1] = d2h[0]/den;
    return jac_inv;
  }

  double operator()(double const& x,double const& y) const{
    return g * h(x,y);
  }

  double max() const{
    unsigned int iteration(0),maxit(1000);
    double error(1.e20),tol(1.e-10),x_const(0.0),y_const(0.0);
    vector<double> x(2,0.e0);
    vector<double> grad(2);
    vector<double> d2h(3);
    vector<vector<double>> jac_inv(2,vector<double>(2));

    // Initialisation des vitesses aux milieu du maillage
    x[0] = 0.5*(b+a); 
    x[1] = 0.5*Ly;

   // Use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){
     // Compute the gradient
     grad = h_prime(x[0],x[1]);
     // Compute the inverse of the jacobian
     jac_inv = h_prime2_inv(x[0],x[1]);
     // Compute new solution
     x[0] = x[0]-jac_inv[0][0]*grad[0] - jac_inv[0][1]*grad[1];
     x[1] = x[1]-jac_inv[1][0]*grad[0] - jac_inv[1][1]*grad[1];
     // Compute new error
     grad = h_prime(x[0],x[1]);
     error = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
     // Increase iterator
     iteration += 1;
   }

   // Send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // Check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;
   } 
   else {
     d2h = h_prime2(x[0],x[1]);  
     if(!(d2h[0]<0 && d2h[0]*d2h[2]-d2h[1]*d2h[1]>0)){
       cout << "Newton method: maximum not found!" << endl;
     }
   }
   // Calculer la vitesse maximale
   x_const = std::max(a,std::min(x[0],b)); 
   y_const = std::max(0.,std::min(x[1],Ly));
   return g * std::max(h0,std::max(h(x_const,y_const),h1));
  };
};



// ============================================================================================================================================================
// Main
int main(int argc, char* argv[]){
  // ----------------------------------------------------------------------------------------------------------------------
  // Déclarations initiales
  unsigned int i,j;
  string inputPath("configuration.in"); // Fichier d'input par défaut
  if(argc>1) // Fichier d'input spécifié par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les paramètres sont lus et stockés dans une "map" de strings.

  for(int k=2; k<argc; ++k) // Input complémentaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[k]);

  // Paramètres de simulation :
  double tfin		 = configFile.get<double>("tfin");
  unsigned int Nx_real	 = configFile.get<int>("Nx");
  unsigned int Ny_real	 = configFile.get<int>("Ny");
  double CFL		 = configFile.get<double>("CFL");
  string type_u2 	 = configFile.get<string>("type_u2").c_str();
  double pert_amplitude = configFile.get<double>("pert_amplitude");
  double pert_velocity  = configFile.get<double>("pert_velocity");
  int mode_num_x   = configFile.get<int>("mode_num_x");
  int mode_num_y	 = configFile.get<int>("mode_num_y");
  bool pertu_f_tilde	 = configFile.get<bool>("pertu_f_tilde");

  // Génération de la fonction u2(x,y)
  U2* u2;
  if(type_u2 == "const")
    u2 = new U2_const(configFile);
  else if(type_u2 == "onde_cas1")
    u2 = new U2_onde_cas1(configFile);
  else if(type_u2 == "onde_cas2")
    u2 = new U2_onde_cas2(configFile);  
  else{
    cerr << "Please select type_u2=""const"" or ""onde_cas1"" or ""onde_cas2""." << endl;
    return -1;
  }

  // ----------------------------------------------------------------------------------------------------------------------
  // TODO : Calculer le maillage et le pas de temps
  double L_x(u2->get_right_extremum()- u2->get_left_extremum());
  double L_y(u2->get_upper_extremum() - u2->get_lower_extremum());

  double dx(L_x/(Nx_real-1));
  double dy(L_y/(Ny_real-1));
  double dt(sqrt(pow(CFL,2)/(u2->max()*((1/pow(dx,2)) + (1/pow(dy,2))))));
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Conditions aux bords (les strings sont converties en valeurs numériques à l'aide d'un énumérateur) :
  bound_cond bc_left, bc_right, bc_upper, bc_lower;
  unsigned int Nx=Nx_real, Ny=Ny_real, start_xiD=0, end_xiD=Nx_real-1, start_yiD=0, end_yiD=Ny_real-1; 

  bc_left = read_boundary_condition(configFile.get<string>("bc_left").c_str());
  bc_right = read_boundary_condition(configFile.get<string>("bc_right").c_str());
  bc_lower = read_boundary_condition(configFile.get<string>("bc_lower").c_str());
  bc_upper = read_boundary_condition(configFile.get<string>("bc_upper").c_str());

  // ----------------------------------------------------------------------------------------------------------------------
  // Lire certaines données si un des bords est harmonique
  double A, omega; // Paramètres d'excitation
  // Lire si une impulsion doit être éxecutée:
  // C'est une option facultative, qui veut dire qu'on excite une seule periode d'excitation
  bool impulsion(configFile.get<bool>("impulsion"));
  bool neumann_apres_impulsion = configFile.get<bool>("neumann_apres_impulsion"); // Impulsion et loi du bord après impulsion
  if(bc_left == harmonic || bc_right == harmonic || bc_upper== harmonic || bc_lower == harmonic){
    A = configFile.get<double>("A");
    omega = configFile.get<double>("omega");
  }

  // ----------------------------------------------------------------------------------------------------------------------
  // Fichiers de sortie :
  bool write_mesh = configFile.get<bool>("write_mesh"); // Exporter x,y
  bool write_f = configFile.get<bool>("write_f"); // Exporter f(x,y,t) ou non

  ofstream fichier_mesh(configFile.get<string>("output_mesh").c_str());
  fichier_mesh.precision(15);
  
  ofstream fichier_f(configFile.get<string>("output_file").c_str());
  fichier_f.precision(15);

  ofstream fichier_E(configFile.get<string>("output_energy").c_str());
  fichier_E.precision(15);

  ofstream fichier_u(configFile.get<string>("output_velocity").c_str());
  fichier_u.precision(15);
  double memoire(0.);
  for(double x(u2->get_left_extremum()); x <= u2->get_right_extremum() + .5*dx; x += dx){
    for(double y(u2->get_lower_extremum()); y <= u2->get_upper_extremum() - .5*dy; y += dy){
      memoire = y;
      fichier_u << sqrt((*u2)(x,y)) << " ";
    }
    fichier_u << sqrt((*u2)(u2->get_right_extremum(), memoire + dy)) << endl;
    memoire = 0;
  }
  fichier_u.close();

  // ----------------------------------------------------------------------------------------------------------------------
  // Lire les données pour initialiser les tableaux
  string type_init(configFile.get<string>("type_init").c_str());
  double F0(configFile.get<double>("F0"));

  // TODO: Calcul du vecteur d'onde selon le mode propre en x et y: k_x=m*pi/L_x; k_y=n*pi/L_y;
  double k_wave_x(mode_num_x*M_PI/L_x);
  double k_wave_y(mode_num_y*M_PI/L_y);

  // Put has first line the position vector
  vector<double> x_mesh(Nx_real),y_mesh(Ny_real);
  for(i = start_xiD; i < end_xiD + 1; ++i){
    x_mesh[i] = u2->get_left_extremum() + (i-start_xiD)*dx; 
  }
  for(i = start_yiD; i < end_yiD + 1; ++i){
    y_mesh[i] = u2->get_lower_extremum() + (i-start_yiD)*dy;
  }
  if(write_mesh){
    fichier_mesh << x_mesh << endl;
    fichier_mesh << y_mesh << endl;
    fichier_mesh.close();
  }
  
  // Initialisation des tableaux du schéma numérique :
  vector<vector<double>> fpast(Nx,vector<double>(Ny)), fnow(Nx,vector<double>(Ny)),\
  fnext(Nx,vector<double>(Ny)), betax2(Nx,vector<double>(Ny)),betay2(Nx,vector<double>(Ny)); 

  if(type_init=="harmonic"){
    // On initialise alors un mode propre (m,n); 
    // TODO: completer fnow, fpast, beta_x^2, beta_y^2
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)
    for(unsigned int i(start_xiD); i <= end_xiD; ++i){
      for(unsigned int j(start_yiD); j <= end_yiD; ++j){
        fpast[i][j] = F0*(cos(x_mesh[i]*k_wave_x - y_mesh[j]*k_wave_y) - cos(x_mesh[i]*k_wave_x + y_mesh[j]*k_wave_y));
        fnow[i][j] = F0*(cos(x_mesh[i]*k_wave_x - y_mesh[j]*k_wave_y) - cos(x_mesh[i]*k_wave_x + y_mesh[j]*k_wave_y));
      }     
    }
  } 
  else{
    // On initialise alors une perturbation nulle.
    // TODO: completer fnow, fpast, beta_x^2, beta_y^2
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)
    for(unsigned int i(start_xiD); i<=end_xiD; ++i){
      for(unsigned int j(start_yiD); j<=end_yiD; ++j){
        fpast[i][j] = 0;
        fnow[i][j] = 0;
      }     
    }
  }

  // ----------------------------------------------------------------------------------------------------------------------
  // Boucle temporelle :
  double t;
  unsigned int stride(0);
  unsigned int n_stride(configFile.get<unsigned int>("n_stride"));

  for(t=0.; t<tfin; t+=dt){
    // Écriture :
    if(stride >= n_stride){
      if(write_f) 
        print_table(fichier_f,t,fnow,start_xiD,start_yiD,end_xiD,end_yiD);
      fichier_E << t << " " << energy(fnow,dx,dy,start_xiD,start_yiD,end_xiD,end_yiD) << endl;
      stride = 0;
    }
    ++stride;

    // TODO: Calculer fnext selon le schéma explicite à 3 niveaux
    for(i=1; i<Nx-1; ++i){
      for(j=1; j<Ny-1; ++j){
        fnext[i][j] = pow(dt,2)*((((*u2)(x_mesh[i+1],y_mesh[j]) - (*u2)(x_mesh[i-1],y_mesh[j]))*(fnow[i+1][j] - fnow[i-1][j])/(4*pow(dx, 2))) \
                    + (((*u2)(x_mesh[i],y_mesh[j+1]) - (*u2)(x_mesh[i],y_mesh[j-1]))*(fnow[i][j+1] - fnow[i][j-1])/(4*pow(dy, 2))) \
                    + (*u2)(x_mesh[i],y_mesh[j])*((fnow[i+1][j] - 2*fnow[i][j] + fnow[i-1][j])/pow(dx, 2))\
                    + (*u2)(x_mesh[i],y_mesh[j])*((fnow[i][j+1] - 2*fnow[i][j] + fnow[i][j-1])/pow(dy, 2)))\
                    + 2*fnow[i][j] - fpast[i][j]; // À modifier!
      } 
    }



    // ----------------------------------------------------------------------------------------------------------------------
    // Conditions aux bords:
    double xL(0); double xR(0); double yL(0); double yU(0);
    switch(bc_left){ // Condition au bord "gauche" (x=0)
      case dirichlet: // TODO : Compléter la condition au bord gauche dirichlet homogène ("fixe")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[start_xiD][y] = fnow[start_xiD][y];
        }
        break;
      case neumann: // TODO : Compléter la condition au bord gauche neumann homogène ("libre")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[start_xiD][y] = fnext[start_xiD+1][y];
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord gauche harmonique 
        if (impulsion && t >= (2*M_PI/omega)){
          for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
            fnext[start_xiD][y] = 0;
          }
          if (neumann_apres_impulsion){
            bc_left = neumann;
          }
          else{
            bc_left = dirichlet;
          }
          break;
        }
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[start_xiD][y] = perturbation(t, xL, y*dy, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, pertu_f_tilde);
        }
        break;
      default:
        throw "Invalid left boundary condition!";
    }

    switch(bc_right){ // Condition au bord droite (x=L_x)
      case dirichlet: // TODO : Compléter la condition au bord droit dirichlet homogène ("fixe")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[end_xiD][y] = fnow[end_xiD][y];
        }
        break;
      case neumann: // TODO : Compléter la condition au bord droit neumann homogène ("libre")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[end_xiD][y] = fnext[end_xiD-1][y];        
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord droit harmonique
        if (impulsion && t >= (2*M_PI/omega)){
          for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
            fnext[end_xiD][y] = 0;
          }
          if (neumann_apres_impulsion){
            bc_right = neumann;
          }
          else{
            bc_right = dirichlet;
          }
          break;
        }
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[end_xiD][y] = perturbation(t, xR, y*dy, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, pertu_f_tilde);
        }
        break;
      default:
        throw "Invalid right boundary condition!";
    }
    
    switch(bc_lower){ // Condition au bord inférieur (y=0)
      case dirichlet: // TODO : Compléter la condition au bord inferieur dirichlet homogène ("fixe")
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][start_yiD] = fnow[x][start_yiD];
        }
        break;
      case neumann:  // TODO : Compléter la condition au bord inferieur neumann homogène ("libre")
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][start_yiD] = fnext[x][start_yiD+1];
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord inferieur harmonique
        if (impulsion && t >= (2*M_PI/omega)){
          for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
            fnext[end_xiD][y] = 0;
          }
          if (neumann_apres_impulsion){
            bc_lower = neumann;
          }
          else{
            bc_lower = dirichlet;
          }
          break;
        }
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][start_yiD] = perturbation(t, x*dx, yL, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, pertu_f_tilde);
        }
        break;
      default:
        throw "Invalid lower boundary condition!";
    }

    switch(bc_upper){ // Condition au bord supérieur (y=L_y)
      case dirichlet: // TODO : Compléter la condition au bord superieur dirichlet homogène ("fixe")
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][end_yiD] = fnow[x][end_yiD];
        }
        break;
      case neumann:  // TODO : Compléter la condition au bord superieur neumann homogène ("libre")
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][end_yiD] = fnext[x][end_yiD-1];
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord superieur harmonique
        if (impulsion && t >= (2*M_PI/omega)){
          for (unsigned int x(start_xiD); x <= end_xiD; ++x){
            fnext[x][end_yiD] = fnow[x][end_yiD];
          }
          if (neumann_apres_impulsion){
            bc_upper = neumann;
          }
          else{
            bc_upper = dirichlet;
          }
          break;
        }
        for (unsigned int x(start_xiD); x <= end_xiD; ++x){
          fnext[x][end_yiD] = perturbation(t, x*dx, yU, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, pertu_f_tilde);
        }
        break;
      default:
        throw "Invalid upper boundary condition!";
    }
    
    // ----------------------------------------------------------------------------------------------------------------------
    // Mise à jour :
    fpast = fnow;
    fnow  = fnext;
  }
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Écrire la solution
  if(write_f) print_table(fichier_f,t,fnow,start_xiD,start_yiD,end_xiD,end_yiD);
    fichier_E << t << " " << energy(fnow,dx,dy,start_xiD,start_yiD,end_xiD,end_yiD) << endl;

  fichier_f.close();
  fichier_E.close();

  return 0;
}
