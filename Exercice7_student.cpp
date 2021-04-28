#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.h"

using namespace std;
typedef enum{not_defined,dirichlet,neumann,harmonic} bound_cond;

// methode pour lire les conditions aux limites
bound_cond read_boundary_condition(string const& bnd_in){
  if(bnd_in == "dirichlet")
    return dirichlet;
  else if(bnd_in == "neumann")
    return neumann;
  else if(bnd_in == "harmonic")
    return harmonic;
  else
  {
    cerr << "Please select "+bnd_in+"=""dirichlet"", ""neumann"", ""harmonic""." << endl;
    return not_defined;
  }
}

//
// Objets pour la fonction u^2(x,y)
//lire les limites au bord
class U2 {
public:

  U2(ConfigFile const& configFile){
    xL=configFile.get<double>("xL"); 	//limite a gauche du fond marin
    xR=configFile.get<double>("xR");	//limite a droite du fond mari
    yL=configFile.get<double>("yL");	//limite inferieur du fond marin
    yU=configFile.get<double>("yU");	//limite superieur du fond marin
  }
//méthodes get
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
  // Methodes virtuelles pures => classe abstraite
  virtual double operator()(double const& x, double const& y) const = 0; // Evalue u^2 au point (x,y)
  virtual double max() const = 0; // Renvoie max(u^2(x,y))

protected: // variables en protected
  double xL,xR,yL,yU;

};

class U2_const: public U2 { //vitesse de propagation au carré (scalaire)
public:
  // Pas de constructeur par defaut => on force a specifier une valeur
  U2_const(ConfigFile const& configFile) :
  U2(configFile), u2(pow(configFile.get<double>("u"),2))		
  {}

  // Definition des methodes virtuelles pures :
  double operator()(double const& x,double const& y) const
  {
    return u2;
  }

  double max() const
  {
    return u2;
  }

private:
  double u2;
};



class U2_onde: public U2 {		//caractérise l'onde
public:

  U2_onde(ConfigFile const& configFile) :
  U2(configFile),
  g(configFile.get<double>("g")),		//acceleration de gravite
  h0(configFile.get<double>("h0")),		//profondeur maximale du fond marin
  h1(configFile.get<double>("h1")),		//profondeur minimale du fond marin
  a(configFile.get<double>("a")),		//limite inferieur de l'onde en x
  b(configFile.get<double>("b"))		//limite superieur de l'onde en y
  {}

protected:
  double g, h0, h1, a, b, Lx;
};

class U2_onde_cas1: public U2_onde {		//equation d'onde dans le 1er cas: Vitesse de propagation constante
public:
  U2_onde_cas1(ConfigFile const& configFile) :
  U2_onde(configFile) {}
  
  double h(double const& x,double const& y) const {  // calcule la profondeur
  double h_val(h0);
  if(x>a && x<b){ // si conditions limites sont respectées
	  h_val = h0+(h1-h0)*sin(M_PI*(x-a)/(b-a)); //M_PI est juste le pi habituel
	  } 
  return h_val;
  }
  
  double h_prime(double const& x) const {
    return (M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a)))/(b-a); //retourne la valeur de h_prime
  }
  
  double h_prime2(double const& x) const {
    return -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a)))/((b-a)*(b-a)); //retourne la valeur de h_2_prime
  }
  
  double operator()(double const& x,double const& y) const
  {
    return g * h(x,y);		// ??????????
  }
  
  double max() const
  {
    unsigned int iteration(0),maxit(1000);  // seulement des int permetant d'avoir la 1er et la dernière itération
    double error(1.e20),tol(1.e-10);		// double pour l'erreur et la tolérance 	l'erreur est très élevée a la première itération pour qu'il y en ait au moins une, elle sera changé juste après
    double x(5.e-1*(a+b));

   // use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){		//on cherche a avoir l'erreur plus basse que la tolérance pour avoir une convergence
     // compute new solution
     x = x - (h_prime(x)/h_prime2(x));
     // compute new error
     error = abs(h_prime(x));
     // increase iterator
     iteration += 1;
   }

   // send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;	//le 1er if implique cela
   } else {
     if(h_prime2(x)<0.e0){
       cout << "Newton method: maximum not found!" << endl;
     }
   }
   // calculer la vitesse maximale
   return g*(std::max(h0,std::max(h(std::max(a,std::min(x,b)),0.),h1)));
  }
  
};



class U2_onde_cas2: public U2_onde {		//equation d'onde dans le 1er cas: Vague avec profondeur de l’oc ́ean variable
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

  vector<double> h_prime(double const& x,double const& y) const
  {
    vector<double> grad(2);		//gradiant
    grad[0] = (M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/(b-a);
    grad[1] = (M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*cos(M_PI*y/Ly))/Ly;
    return grad;
  }

  vector<double> h_prime2(double const& x,double const& y) const {
    vector<double> d2h(3,0.0);
    d2h[0] = -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/((b-a)*(b-a));
    d2h[1] = (M_PI*M_PI*(h1-h0)*cos(M_PI*(x-a)/(b-a))*cos(M_PI*y/Ly))/((b-a)*Ly);
    d2h[2] = -(M_PI*M_PI*(h1-h0)*sin(M_PI*(x-a)/(b-a))*sin(M_PI*y/Ly))/(Ly*Ly);
    return d2h;
  }

  vector<vector<double>> h_prime2_inv(double const& x,double const& y) const
  {
    vector<vector<double>> jac_inv(2,vector<double>(2));
    vector<double> d2h(h_prime2(x,y));
    double den=0.0;
    
    // calculer l inverse du jacobian
    den = d2h[0]*d2h[2]-d2h[1]*d2h[1];
    jac_inv[0][0] = d2h[2]/den;
    jac_inv[0][1] = -d2h[1]/den;
    jac_inv[1][0] = d2h[1];
    jac_inv[1][1] = d2h[0]/den;
    return jac_inv;
  }

  double operator()(double const& x,double const& y) const
  {
    return g * h(x,y);
  }

  double max() const
  {

    unsigned int iteration(0),maxit(1000);
    double error(1.e20),tol(1.e-10),x_const(0.0),y_const(0.0);
    vector<double> x(2,0.e0);
    vector<double> grad(2);
    vector<double> d2h(3);
    vector<vector<double>> jac_inv(2,vector<double>(2));

    // initialisation des vitesses aux milieu du maillage
    x[0] = 0.5*(b+a); 
    x[1] = 0.5*Ly;

   // use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){			//on teste d'aitres équations que dans le 1er cas
     // compute the gradient
     grad = h_prime(x[0],x[1]);
     // compute the inverse of the jacobian
     jac_inv = h_prime2_inv(x[0],x[1]);
     // compute new solution
     x[0] = x[0]-jac_inv[0][0]*grad[0] - jac_inv[0][1]*grad[1];
     x[1] = x[1]-jac_inv[1][0]*grad[0] - jac_inv[1][1]*grad[1];
     // compute new error
     grad = h_prime(x[0],x[1]);
     error = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
     // increase iterator
     iteration += 1;
   }

   // send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;
   } else {
       d2h = h_prime2(x[0],x[1]);  
       if(!(d2h[0]<0 && d2h[0]*d2h[2]-d2h[1]*d2h[1]>0)){
         cout << "Newton method: maximum not found!" << endl;
       }
     }
   // calculer la vitesse maximale
   x_const = std::max(a,std::min(x[0],b)); 
   y_const = std::max(0.,std::min(x[1],Ly));
   return g * std::max(h0,std::max(h(x_const,y_const),h1));
  };
};

//
// TODO : Calcul de l'energie de l'onde
//

double energy(vector<vector<double>> const& f, double const& dx, double const& dy,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD)
{
  double energy_(0.);
  // TODO: completer ici en utilisant la formule des trapezes
  double sum_f(0.0);
  for(double i(start_xiD); i<end_xiD ; i+=1){ //il a pas besoin d'être int et commencer a 0
	  for(double j(start_yiD); j< end_yiD ; j+=1){
	  sum_f = (f[i][j] * f[i][j] + f[i+1][j] * f[i+1][j] + f[i][j+1] * f[i][j+1] + f[i+1][j+1] * f[i+1][j+1]);  //équation donnant l'énergie totale
	}
  }
energy_ = (1/4) * dx * dy * sum_f;

  return energy_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau 1D
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)		
{
  unsigned int len(v.size());

  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";

  if(len > 0)
    o << v[len-1];

  return o;
}
//
// Surcharge de l'operateur pour ecrire les elements d'un tableau 2D
//
template<typename T> void print_table(ostream& o, double const& t , vector<vector<T>> const& v,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD)
{
  unsigned int i,j;
  
  for(i=start_yiD; i<end_yiD+1; ++i){
      o << t << " ";
      for(j=start_xiD; j<end_xiD; ++j){
        o << v[i][j] << " "; 
      }
      o << v[i][end_xiD] << endl;
  }
}

// TODO: calculer la perturbation
double perturbation(double const& t,double const& x,double const& y,\
double const& xL, double const& xR,double const& yL, double const& yU, \
double const& pert_velocity, double const& pert_amplitude, \
int const mode_num_x,int const mode_num_y, bool const& simple) {
 if (simple){
    return pert_amplitude*cos(M_PI*((mode_num_x*x/(xR - xL))) + (mode_num_y*y/(yU - yL))); // Acos(pi x/L_x + pi y/L_y)
  }
  else{
    return pert_amplitude*sin(pert_velocity*t); //a_e sin(wt)
  }

 return 0.e0;

}

//
// Main
//
int main(int argc, char* argv[])
{
  unsigned int i,j;
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int k=2; k<argc; ++k) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[k]);

  // Parametres de simulation :
  double tfin		 = configFile.get<double>("tfin");
  unsigned int Nx_real	 = configFile.get<int>("Nx");
  unsigned int Ny_real	 = configFile.get<int>("Ny");
  double CFL		 = configFile.get<double>("CFL");
  string type_u2 	 = configFile.get<string>("type_u2").c_str();
  double pert_amplitude = configFile.get<double>("pert_amplitude");
  double pert_velocity  = configFile.get<double>("pert_velocity");
  int mode_num_x	 = configFile.get<int>("mode_num_x");
  int mode_num_y	 = configFile.get<int>("mode_num_y");
  bool simple	 = configFile.get<bool>("simple");

  U2* u2;
  if(type_u2 == "const")
    u2 = new U2_const(configFile);
  else if(type_u2 == "onde_cas1")
    u2 = new U2_onde_cas1(configFile);
  else if(type_u2 == "onde_cas2")
    u2 = new U2_onde_cas2(configFile);  
  else
  {
    cerr << "Please select type_u2=""const"" or ""onde_cas1"" or ""onde_cas2""." << endl;
    return -1;
  }

  // TODO calculer le maillage et le pas de temps

  bool write_mesh = configFile.get<bool>("write_mesh"); // Exporter x,y
  bool write_f = configFile.get<bool>("write_f"); // Exporter f(x,y,t) ou non
  double L_x(u2->get_right_extremum()- u2->get_left_extremum()); //u2 --> pointeur
  double L_y(u2->get_upper_extremum() - u2->get_lower_extremum());

  double dx(L_x/(Nx_real-1));
  double dy(L_y/(Ny_real-1));
  double dt(sqrt(pow(CFL,2)/(u2->max()*((1/pow(dx,2)) + (pow(dy,2))))));


  //Conditions aux bords (les strings sont converties en valeurs numeriques a l'aide d'un enumerateur) :
  bound_cond bc_left, bc_right, bc_upper, bc_lower;
  unsigned int Nx=Nx_real,Ny=Ny_real,start_xiD=0,end_xiD=Nx_real-1,start_yiD=0,end_yiD=Ny_real-1;
  bc_left = read_boundary_condition(configFile.get<string>("bc_left").c_str());		
  if(bc_left==neumann){
    // TODO
  }
  bc_right = read_boundary_condition(configFile.get<string>("bc_right").c_str());
  if(bc_right==neumann) ; // TODO
  bc_lower = read_boundary_condition(configFile.get<string>("bc_lower").c_str());
  if(bc_lower==neumann){
    // TODO
  }
  bc_upper = read_boundary_condition(configFile.get<string>("bc_upper").c_str());
  if(bc_upper==neumann) ; // TODO

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
  
  
  
  // Fichiers de sortie :
  ofstream fichier_mesh(configFile.get<string>("output_mesh").c_str());		//nome du fichier d'output du maillage
  fichier_mesh.precision(15);
  
  ofstream fichier_f(configFile.get<string>("output_file").c_str());		//nome du fichier d'output de la solution
  fichier_f.precision(15);

  ofstream fichier_E(configFile.get<string>("output_energy").c_str());		//nome du fichier d'output de l'energie
  fichier_E.precision(15);

  ofstream fichier_u(configFile.get<string>("output_velocity").c_str());		//nome du fichier d'output de la vitesse
  fichier_u.precision(15);
  for(double y(u2->get_lower_extremum()); y<=u2->get_upper_extremum()+.5*dy; y+=dy){
    for(double x(u2->get_left_extremum()); x<=u2->get_right_extremum()-.5*dx; x+=dx){
      fichier_u << sqrt((*u2)(x,y)) << " ";
    }
    fichier_u << sqrt((*u2)(u2->get_right_extremum(),y)) << endl;
  }
  fichier_u.close();


  // lire les donnees pour initialiser les tableaux
  string type_init(configFile.get<string>("type_init").c_str());
  double F0(configFile.get<double>("F0"));
  //double u2_loc(0);

  // TODO: calcul de l'inverse de la longueur d'onde en x et y: k_x=m*pi/L_x; k_y=n*pi/L_y;
  double k_wave_x(0.);
  double k_wave_y(0.);

  // Initialisation des tableaux du schema numerique :
  vector<vector<double>> fpast(Ny,vector<double>(Nx)), fnow(Ny,vector<double>(Nx)),\
  fnext(Ny,vector<double>(Nx)), betax2(Ny,vector<double>(Nx)),betay2(Ny,vector<double>(Nx)); 

  if(type_init=="harmonic"){
    // On initialise alors un mode propre (m,n); 
    // TODO: completer fnow, fpast, beta_x^2, beta_y^2
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)
    
    //beta = u delta t / delta x
    
    for(unsigned int i(start_xiD); i<=end_xiD; ++i){
      for(unsigned int j(start_yiD); j<=end_yiD; ++j){
        fpast[i][j] = F0*cos(i*dx*k_wave_x + j*dy*k_wave_y); //F0 = A , f^ = Re(e^(ikx))= cos (1/pi(x*k_x + y*k_y)) , k_x = m*pi/L_x , n*pi/L_y  , x = i*dx , j = j*dy
        fnow[i][j] = F0*cos(i*dx*k_wave_x + j*dy*k_wave_y);  
      }     
    }
    
  } else {
    // On initialise alors une perturbation nulle.
    // TODO: completer fnow, fpast
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)

	for(unsigned int i(start_xiD); i<=end_xiD; ++i){
      for(unsigned int j(start_yiD); j<=end_yiD; ++j){
        fpast[i][j] = 0;
        fnow[i][j] = 0;
      }     
    }


  }

  // Boucle temporelle :
  double t;
  unsigned int stride(0);
  unsigned int n_stride(configFile.get<unsigned int>("n_stride"));
  // put has first line the position vector
  vector<double> x_mesh(Nx_real),y_mesh(Ny_real);
  for(i=start_xiD; i<end_xiD+1; ++i){
    x_mesh[i] = u2->get_left_extremum()+(i-start_xiD)*dx; 
  }
  for(i=start_yiD; i<end_yiD+1; ++i){
    y_mesh[i] = u2->get_lower_extremum()+(i-start_yiD)*dy;
  }
  if(write_mesh){
  
    fichier_mesh << x_mesh << endl;
    fichier_mesh << y_mesh << endl;
    fichier_mesh.close();
  }

  for(t=0.; t<tfin; t+=dt)
  {
    // Ecriture :
    if(stride >= n_stride )
    {
      if(write_f) print_table(fichier_f,t,fnow,start_xiD,start_yiD,end_xiD,end_yiD);
      fichier_E << t << " " << energy(fnow,dx,dy,start_xiD,start_yiD,end_xiD,end_yiD) << endl;
      stride = 0;
    }
    ++stride;
    // TODO: calculer fnext selon le schema explicite a 3 niveaux
       for(i=1; i<Nx-1; ++i){
      for(j=1; j<Ny-1; ++j){
        fnext[i][j] = (dt*dt)*((((*u2)(x_mesh[i+1],y_mesh[j]) - (*u2)(x_mesh[i-1],y_mesh[j]))*(fnow[i+1][j] - fnow[i-1][j])/(4*pow(dx, 2))) 
                    + (((*u2)(x_mesh[i],y_mesh[j+1]) - (*u2)(x_mesh[i],y_mesh[j-1]))*(fnow[i][j+1] - fnow[i][j-1])/(4*pow(dy, 2))) 
                    + (*u2)(x_mesh[i],y_mesh[j])*((fnow[i+1][j] - 2*fnow[i][j] + fnow[i-1][j])/pow(dx, 2))  
                    + (*u2)(x_mesh[i],y_mesh[j])*((fnow[i][j+1] - 2*fnow[i][j] + fnow[i][j-1])/pow(dy, 2)))
                    + 2*fnow[i][j] 
                    - fpast[i][j]; // À modifier!
      } 
    }



/*/

\begin{equation*}
    f(x_i,y_j,t_{k+1}) = h_{t}^2 \bigg[ \underbrace{ u^2(x_i,y_j) \ \frac{f(x_{i+1},y_j) - 2f(x_{i},y_{j}) + f(x_{i-1},y_j)}{h_x^2}}_{u^2 \frac{\partial^2 f}{\partial x^2}}
\end{equation*}
\begin{equation*}
      + \underbrace{u^2(x_i,y_j) \  \frac{f(x_{i},y_{j+1}) - 2f(x_{i},y_{j}) + f(x_{i},y_{j-1})}{h_y^2}}_{u^2 \frac{\partial ^2 f}{\partial y^2}}    +  \underbrace{\bigg( \frac{(u^2(x_{i+1},y_j) - u^2(x_{i-1},y_j)) (f(x_{i+1},y_j) - f(x_{i-1},y_j))}{4h_x^2} \bigg)}_{\frac{\partial u^2}{\partial x} \frac{\partial f}{\partial x}}
\end{equation*}
\begin{equation*}
      + \underbrace{\bigg( \frac{(u^2(x_{i},y_{j+1}) - u^2(x_{i},y_{j-1})) (f(x_{i},y_{j+1}) - f(x_{i},y_{j-1}))}{4h_y^2} \bigg)}_{\frac{\partial u^2}{\partial y} \frac{\partial f}{\partial y}} + a_{e}(x,y)\sin(\omega t) \bigg] + 2f(x_i,y_j,t_{k}) - f(x_i,y_j,t_{k-1})
\end{equation*}
/*/



   // Conditions aux bords:
    double xL(0); double xR(0); double yL(0); double yU(0);
    switch(bc_left){ // Condition au bord "gauche" (x=0)
      case dirichlet: // TODO : Compléter la condition au bord gauche dirichlet homogène ("fixe")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[start_xiD][y] = fnow[start_xiD][y]; //f(x0,y)
        }
        break;
      case neumann: // TODO : Compléter la condition au bord gauche neumann homogène ("libre")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[start_xiD][y] = fnext[start_xiD+1][y]; // f(x0,y) = f(x1,y)
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord gauche harmonique 
        if (impulsion && t >= (2*M_PI/omega)){ //T
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
          fnext[start_xiD][y] = perturbation(t, 0, y*dy, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, simple); //perturbation
        }
        break;
      default:
        throw "Invalid left boundary condition!";
    }

    switch(bc_right){ // Condition au bord droite (x=L_x)
      case dirichlet: // TODO : Compléter la condition au bord droit dirichlet homogène ("fixe")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[end_xiD][y] = fnow[end_xiD][y]; //f(Lx,y,t) = f(Lx,y,t+1)
        }
        break;
      case neumann: // TODO : Compléter la condition au bord droit neumann homogène ("libre")
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
          fnext[end_xiD][y] = fnext[end_xiD-1][y]; //f(x_Lx, y,t+1) = f(x_Lx-1 , y,t+1))        
        }
        break;
      case harmonic: // TODO : Compléter la condition au bord droit harmonique
        if (impulsion && t >= (2*M_PI/omega)){
          for (unsigned int y(start_yiD+1); y < end_yiD; ++y){
            fnext[end_xiD][y] = 0; //f(x_L,y) = 0
          }
          if (neumann_apres_impulsion){ //on reviens au stade normal
            bc_right = neumann;
          }
          else{
            bc_right = dirichlet;
          }
          break;
        }
        for (unsigned int y(start_yiD+1); y < end_yiD; ++y){ //harmonique: impulsion
          fnext[end_xiD][y] = perturbation(t, L_x, y*dy, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y,simple);
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
          fnext[x][start_yiD] = perturbation(t, x*dx, 0, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y, simple);
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
          fnext[x][end_yiD] = perturbation(t, x*dx, L_y, xL, xR, yL, yU, omega, A, mode_num_x, mode_num_y,simple);
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
