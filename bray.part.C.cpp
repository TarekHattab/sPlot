// To enable the functionality provided by Armadillo's various macros,

// simply include them before you include the RcppArmadillo headers.

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

#include <RcppArmadilloExtensions/sample.h>

#include <bigmemory/BigMatrix.h>

#include <Rcpp.h> 

using namespace Rcpp;

using namespace arma;


void BrayPart(const arma::Mat<double>& inBigMat, arma::Mat<double> outBigMat) {
  
  int nRows = inBigMat.n_rows;
  
  int nCols = inBigMat.n_cols;
  
  cout <<  nRows << std::endl;
  
  cout <<  nCols << std::endl;
	
	for(int i = 0; i < nRows; i++){
		
		for(int j = 0; j < nRows; j++){
			
			NumericVector pMin(nCols,0.0);
			
			double sumPi = 0.0 ;
			
			double sumPj = 0.0 ;
			  
			  for (int k = 0; k < nCols; k++){
					
					sumPi += inBigMat(i,k);
					
					sumPj += inBigMat(j,k);
					
					if (inBigMat(i,k) < inBigMat(j,k)) {
					
					pMin[k] = inBigMat(i,k);
					
					}
					
					else {
					
					pMin[k] = inBigMat(j,k);	
					
					}
				
				} 
			
			double aComp = Rcpp::sum(pMin);
			
			double bComp = sumPi - aComp;
			
			double cComp = sumPj - aComp;
			
			double minBC = 0.0;
			 
			 if (bComp < cComp) { 
			 
			 minBC = bComp;
			 
			 }
			 
			 else {
			 
			 minBC = cComp; 
			 
			 }
			 
			 outBigMat(i,j) = minBC / (aComp + minBC);
			 
			 //outBigMat(j,i) = ( bComp  + cComp) / (2 * aComp + bComp + cComp);	
		
		}
	
	}
	
	outBigMat =   outBigMat.t();

}


// [[Rcpp::export]]

void BigBray(SEXP pInBigMat, SEXP pOutBigMat) {

	// First we tell Rcpp that the object we've been given is an external pointer.
 
	XPtr<BigMatrix> xpMat(pInBigMat);

	XPtr<BigMatrix> xpOutMat(pOutBigMat);
  
	BrayPart(arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),

			 arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
	);
}

