#include "kadath_polar.hpp"
#include "mpi.h"
#include "magma_interface.hpp"


using namespace Kadath ;

int main(int argc, char** argv) {

	int rc = MPI_Init(&argc, &argv) ;
	int rank = 0 ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
    
#ifdef ENABLE_GPU_USE
    if(rank==0)
    {
        TESTING_CHECK(magma_init());
        magma_print_environment();
    }
#endif

	int kk ;
	double omega ;
  
	FILE* fin = fopen (argv[1], "r") ;
	Space_polar space(fin) ;
	fread_be (&kk, sizeof(int), 1, fin) ;
	fread_be (&omega, sizeof(double), 1, fin) ;
	Scalar nu (space,fin) ;
	Scalar incA (space, fin) ;
	Scalar incB (space, fin) ;
	Scalar incbt (space, fin) ;
	Scalar phi (space, fin) ;
	fclose(fin) ;

	
	int ndom = space.get_nbr_domains() ;
	double qpi = 4*M_PI ;


	Param_tensor parameters ;
	parameters.set_m_quant() = kk ;
        phi.set_parameters() = parameters ;

	int number = 8 ;
	double step = 0.003 ; 
	for (int kant=0 ; kant<number ; kant++) {
		if (kant !=0)
		omega -= step ;

	if (rank==0)
	  cout << "Computation with omega = " << omega << endl ;

      Scalar one (space) ;
      one = 1 ;
      one.std_base() ;
      Scalar rsint (space) ;
    for (int d=0 ; d<ndom-1 ; d++)
		rsint.set_domain(d) = space.get_domain(d)->mult_r(space.get_domain(d)->mult_sin_theta(one(d))) ;
    rsint.set_domain(ndom-1) = space.get_domain(ndom-1)->mult_sin_theta(one(ndom-1)) ;
  

      System_of_eqs syst (space, 0, ndom-1) ;
  
      syst.add_var ("incB", incB) ;
      syst.add_var ("nu", nu) ;
      syst.add_var ("incA", incA) ;
      syst.add_var ("phi", phi) ;
      syst.add_var ("incbt",incbt) ;
      syst.add_cst ("ome", omega) ;

      syst.add_cst ("rsint", rsint) ;
      syst.add_cst ("k", kk) ;
      syst.add_cst ("qpi", qpi) ;

      syst.add_def ("phisurrsint = divrsint(phi)") ;
      syst.add_def ("ap = exp(nu)") ;
      syst.add_def ("bt = divrsint(incbt)") ;
      syst.add_def ("B = (divrsint(incB) + 1)/ap") ;
      syst.add_def ("A = exp(incA - nu)") ;
      
      syst.add_def ("E = 0.5*(ome-bt*k)^2/ap^2*phi^2 + 0.5*scal(grad(phi),grad(phi))/A^2 + 0.5*phi^2 + 0.5*k*k*phisurrsint*phisurrsint/B^2") ;
      syst.add_def ("Pp = k/ap* (ome-bt*k)*phi^2") ;
      syst.add_def ("S = -0.5*scal(grad(phi), grad(phi))/A^2 - 0.5*k*k*phisurrsint*phisurrsint/B^2 + 1.5*(ome-bt*k)^2/ap^2 * phi^2- 1.5*phi^2") ;
      syst.add_def ("Spp = 0.5*(ome-bt*k)^2/ap^2*phi^2 -0.5* scal(grad(phi),grad(phi))/A^2 -0.5* phi^2 +0.5* k*k*phisurrsint*phisurrsint/B^2") ;
       
	  for (int d=0 ; d<ndom-1 ; d++)
	  syst.add_def (d, "eqB = lap2(incB) -2*qpi*ap*A^2*B*rsint*(S-Spp)") ;
      syst.add_def (ndom-1, "eqB = lap2(incB) -2*qpi*ap*A^2*B*rsint*multr(S-Spp)") ;



	for (int d=0 ; d<ndom-1 ; d++)
	syst.add_def (d, "eqshift = lap(incbt) - divrsint(bt) - rsint*scal(grad(bt),grad(nu-3*log(B)))+4*qpi*ap*A^2/B^2*divrsint(Pp)") ;
      syst.add_def (ndom-1, "eqshift = lap(incbt) - divrsint(bt) - rsint*scal(multr(grad(bt)),grad(nu-3*log(B)))+4*qpi*ap*A^2/B^2*divrsint(Pp)") ;
     
 
      for (int d=0 ; d<ndom-1 ; d++)
	syst.add_def (d, "eqnu = lap(nu) - B^2*rsint^2/2/ap^2 * scal(grad(bt) , grad(bt)) + scal(grad(nu), grad(nu+log(B))) - qpi*A^2*(E+S)") ;
      syst.add_def (ndom-1, "eqnu = lap(nu) - B^2*rsint^2/2/ap^2 * scal(multr(grad(bt)) , multr(grad(bt))) + scal(grad(nu), grad(nu+log(B))) - qpi*A^2*(E+S)") ;
	
    
     for (int d=0 ; d<ndom-1 ; d++) 
	syst.add_def (d, "eqA = lap2(incA) - 2*qpi*A^2*Spp- 3 * B^2 * rsint^2 /4 / ap^2 * scal(grad(bt) , grad(bt)) + scal(grad(nu),grad(nu))") ;
       syst.add_def (ndom-1, "eqA =lap2(incA) - 2*qpi*A^2*Spp- 3 * B^2 * rsint^2 /4 / ap^2 * scal(multr(grad(bt)) , multr(grad(bt))) + scal(grad(nu),grad(nu)) ") ;
    
 
	syst.add_def ("print = lap2(incA)") ;

          for (int d=0 ; d<ndom-1 ; d++)
	syst.add_def (d, "eqphi = lap(phi) - A^2*(1-ome*ome/ap^2+2*bt/ap^2*ome*k-bt^2/ap^2*k*k)*phi + scal(grad(phi),grad(nu+log(B))) - divrsint(A^2/B^2 -1)*divrsint(phi)*k*k") ;
      syst.add_def (ndom-1, "eqphi = lap(phi) - A^2*(1-ome*ome/ap^2+2*bt/ap^2*ome*k-bt^2/ap^2*k*k)*phi + scal(grad(phi),grad(nu+log(B))) - k*k*divrsint(A^2/B^2-1)*divrsint(phi)") ;
  
      
      space.add_eq (syst, "eqB=0", "incB", "dn(incB)") ; 
      space.add_eq (syst, "eqA=0", "incA", "dn(incA)") ; 
      space.add_eq (syst, "eqphi=0", "phi", "dn(phi)") ; 
      space.add_eq (syst, "eqnu=0", "nu", "dn(nu)") ; 
      space.add_eq (syst, "eqshift=0", "incbt", "dn(incbt)") ; 
       
    
      syst.add_eq_bc (ndom-1, OUTER_BC, "nu=0") ;
      syst.add_eq_bc (ndom-1, OUTER_BC, "incA=0") ;
      syst.add_eq_bc (ndom-1, OUTER_BC, "incB=0") ;
      syst.add_eq_bc (ndom-1, OUTER_BC, "incbt=0") ;
      syst.add_eq_bc (ndom-1, OUTER_BC, "phi=0") ;

      double conv ;
      bool endloop = false ;
      int ite = 1 ;
      while (!endloop) {
	endloop = syst.do_newton(1e-8, conv) ;
	if(rank==0)
	        cout << "Newton iteration " << ite << " " << conv  << endl ;
	ite++ ;
	 }


	
	if (rank==0) {
	
		char name[100] ;
		sprintf (name, "bos_%d_%f.dat", kk, omega) ;
		FILE* fiche = fopen (name, "w") ;
		space.save(fiche) ;
		fwrite_be (&kk, sizeof(int), 1, fiche) ;
		fwrite_be (&omega, sizeof(double), 1, fiche) ;
		nu.save(fiche) ;
		incA.save(fiche) ;
		incB.save(fiche) ;
		incbt.save(fiche) ;
		phi.save(fiche) ;
		fclose(fiche) ;
		}
	}

#ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_finalize());
	}
#endif
    MPI_Finalize() ;
	
}
