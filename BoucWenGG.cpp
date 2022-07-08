/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 0.1
// $Date: 2021/02/19
// $Source: 
                                                                        
// Written: Andrea Lucchini (andrea.lucchini@uniroma1.it)
//
// Description: This file contains the class implementation for the
// BoucWenGG Uniaxial Material (refer to Gerolymos and Gazetas, 2005, and Drosos et al., 2012).
//
// Modified: Andrea Marchi (andrea.marchi@uniroma1.it) [19/02/2021]
// 

#include <elementAPI.h>
#include "BoucWenGG.h"
#include <MaterialResponse.h>
#include <Information.h>



// External procedure to parse the material command line
#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numBoucWenGG = 0;

OPS_Export void* OPS_BoucWenGG()
{

  // print out some KUDO's
  if (numBoucWenGG == 0) {
    opserr << "BoucWenGG unaxial material (Drosos et al., 2012) - Written by Andrea Lucchini,\nmodified by Andrea Marchi" << endln;
    numBoucWenGG =1;
  }

  // pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // parse the input line for the material parameters
  int    iData[2];
  double dData[8];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING error reading uniaxialMaterial BoucWenGG tag" << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading uniaxialMaterial BoucWenGG data (alpha, ko, epsY, n, gamma, beta, s1, s2)\nThey have to be double values!" << endln;
    return 0;	
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING error reading uniaxialMaterial BoucWenGG data (mkur)\nIt has to be an integer value!" << endln;
    return 0;
  }
  
  numData = 1;
  double tol;
  if (OPS_GetDoubleInput(&numData, &tol) != 0) {
    tol = 1.0e-8;
  }
  
  numData = 1;
  int steps;
  if (OPS_GetIntInput(&numData, &steps) != 0) {
    steps = 20;
  }

   // create a new material
  theMaterial = new BoucWenGG(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], iData[1], tol, steps);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BoucWenGG" << endln;
    return 0;
  }

  // return the material
  return theMaterial;

}



BoucWenGG::BoucWenGG(int tag, 
					double p_alpha,
					double p_ko,
					double p_strainY,
					double p_n,
					double p_gamma,
					double p_beta,
					double p_s1,
					double p_s2,
					int p_mkur,
					double p_tolerance,
					int p_maxNumIter)
:UniaxialMaterial(tag,0),
alpha(p_alpha), ko(p_ko), strainY(p_strainY), n(p_n), gamma(p_gamma), beta(p_beta), s1(p_s1), s2(p_s2), mkur(p_mkur),
tolerance(p_tolerance), maxNumIter(p_maxNumIter)
{
	// Initialize variables
    this->revertToStart();
}



BoucWenGG::~BoucWenGG()
{
	// Nothing to do here
}



double BoucWenGG::signum(double value)
{
	//sign function [sgn(x)]
	if (value > 0.0)       return 1.0;
	else if (value < 0.0)  return -1.0;
	else                   return 0.0; //if x=0 sgn(x)=0 [avoid computational errors]
}




double BoucWenGG::fun(double z, double z_old) {
	//are all global variables for the costitutive model (gamma, beta, n,  TdStrain, Ttheta)
	return z - z_old - (1.0 - pow(fabs(z), n) * (gamma + beta * (signum(TdStrain * z)))) * Ttheta * TdStrain / strainY;
}



int BoucWenGG::setTrialStrain (double strain, double strainRate)
{
	
	// Set trial strain and compute trial strain increment
	Tstrain = strain;
	TdStrain = Tstrain - Cstrain;
	
	// Initial declarations (make sure not to declare class variables here!)
	double Psi, Phi, f, Phi_, f_, Tznew, Tzold, sign;
	double dStress;	// stress increment
	
	
	// Compute deterioration parameter (theta)
	// a) check reversal:
	// - in the positive stress domain
	if ( (signum(CdStrain*TdStrain) < 0.0 ) && (Cstrain == CmaxStrain ) ) {
		TNrevp = CNrevp + 1;									// update counter of reversals
		TmaxRevStrain = Cstrain;								// update value of the attained max strain at reversal
		if ( TNrevp == 1 ) Tmurp = 0.0;							// at the first reversal the deterioration parameter is equal to 0
		else Tmurp = (TmaxRevStrain-CmaxRevStrain)/(2*strainY);
	} else {
		TmaxRevStrain = CmaxRevStrain;
		Tmurp = Cmurp;
	}
	// - in the negative stress domain
	if ( (signum(CdStrain*TdStrain) < 0.0 ) && (Cstrain == CminStrain ) ) {
		TNrevn = CNrevn + 1;									// update counter of reversals
		TminRevStrain = Cstrain;								// update value of the attained min strain at reversal
		if ( TNrevn == 1 ) Tmurn = 0.0;							// at the first reversal the deterioration parameter is equal to 0
		else Tmurn = (CminRevStrain-TminRevStrain)/(2*strainY);
	} else {
		TminRevStrain = CminRevStrain;
		Tmurn = Cmurn;
	}
	// b) calculate the reference strain ductility
	if ( Tmurp > Tmurn ) Tmur = Tmurp; else Tmur = Tmurn;
	// c) calculate the deterioration parameter
	if ( Tmur < s2 ) Ttheta = 1.0; else Ttheta = (s1+alpha*(Tmur-1)+s2)/(s1+Tmur);
	
	
	// Newton-Raphson scheme to solve for z_{i+1} := Tznew
	int count = 0;
	double startPoint = 00; //prima era 0.01 poi avevo messo 0.001 (mi sa che non può essere zero per il calcolo della derivata...)    METTERE Cz ????????
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;
	while ( (fabs(Tzold-Tznew) > tolerance) && count<maxNumIter) {

		sign = signum(TdStrain*Tz);
		Psi = gamma + beta*sign;
		Phi = 1.0 - pow(fabs(Tz),n)*Psi;
		f = Tz - Cz - Phi*Ttheta*TdStrain/strainY;

		// Evaluate function derivative f' (underscore:=prime)
		sign = signum(Tz);
		double pow1;
		if (Tz == 0.0) pow1 = 0.0; else pow1 = pow(fabs(Tz),(n-1));
		Phi_ = - n*pow1*sign*Psi;
		f_ = 1.0 - Phi_*Ttheta*TdStrain/strainY;

		// Issue warning if derivative is zero
		if ( fabs(f_)<1.0e-10 ) opserr << "WARNING: BoucWenGG::setTrialStrain() -- zero derivative in Newton-Raphson scheme" << endln;

		// Take a Newton step
		Tznew = Tz - f/f_;
		//if (fabs(Tznew) > 2 * fabs(Tz) && Tz != 0) { Tznew = 2 * fabs(Tz) * signum(Tznew); opserr << "Limitazione eseguita!!!!" << endln; } //limita il "salto" che può fare...
		//if (signum(Tznew) != signum(Tz)) opserr << "Inversione di segno!!" << endln;

		// Update the root (but keep the old value for convergence check)
		Tzold = Tz; 
		Tz = Tznew;

		// Update counter
		count++;

		// If Newton-Raphson take too long it means that is not converging... (MIGLIORARE!)
		if (count == 50) {
			// take some bisection steps:
			//opserr << "bisezione... ("<< fabs(Tzold - Tznew) <<") ";
			int contatore = 0;
			// search for correct interval
			double a = 0, b; //interval extrema
			double max_z = 2*pow(1 / (gamma + beta), 1/n); //maximum value for z
			if (fun(a, Cz) * fun(max_z, Cz) < 0) { b = max_z; /*opserr << "(pos)  ";*/ }
			else { b = -max_z; /*opserr << "(neg)  ";*/ }
			while ((fabs(a - b) > tolerance*10) && (contatore < 100)) {
				double c = 0.5 * (a + b); //interval midpoint
				if (fun(a, Cz) * fun(c, Cz) < 0) b = c; //correct interval is [a,c]
				else a = c; //correct interval is [c,b]
				contatore++;
			}
			Tz = 0.5 * (a + b);
			//opserr << "found Tz = " << Tz << "\t [contatore = " << contatore << ", count = " << count << "]";
			if (fabs(a - b) > tolerance) { count++; /*opserr << " forcing steps\n";*/ } //forcing other steps
			else { Tznew = Tzold = Tz; /*opserr << " forcing exit\n";*/ } //have reached convergence, forcing loop exit
		}

		// If steps are too much algorithm try to finish with "false position" algorithm
		if (count == 100) {
			// find solution with "false position" algorithm:
			//opserr << "falsa posizione... ("<< fabs(Tzold - Tznew) <<") ";
			int contatore = 0;
			// search for correct interval
			double a = 0, b; //interval extrema
			double max_z = 2*pow(1/(gamma+beta),1/n); //maximum value for z
			if (fun(a, Cz) * fun(max_z, Cz) < 0) { b = max_z; /*opserr << "(pos)  ";*/ }
			else { b = -max_z; /*opserr << "(neg)  ";*/ }
			while ((fabs(a - b) > tolerance) && (contatore < 100)) {
				double c = a - fun(a, Cz) * (b - a) / (fun(b, Cz) - fun(a, Cz)); //interval false position
				if (fun(a, Cz) * fun(c, Cz) < 0) b = c; /*correct interval is [a,c]*/ else a = c; /*correct interval is [c,b]*/
				contatore++;
			}
			Tz = a - fun(a, Cz) * (b - a) / (fun(b, Cz) - fun(a, Cz));
			//opserr << "found Tz = " << Tz << "\t [contatore = " << contatore << ", count = " << count << "]";
			Tznew = Tzold = Tz; //have reached convergence, forcing loop exit
		}

		// Issue warning if we didn't converge
		//opserr << count << ":  fabs(Tzold-Tznew) = fabs(" << Tzold << "-" << Tznew << ") = " << fabs(Tzold - Tznew) << endln;
		if (count >= maxNumIter) {
			opserr << "WARNING: BoucWenGG::setTrialStrain() -- did not find the root z_{i+1}, after " << maxNumIter << " iterations and norm: " << fabs(Tzold-Tznew) << " (max: " << tolerance << ")" << endln;
		}

		// Compute stress
		Tstress = alpha*ko*Tstrain + (1-alpha)*ko*strainY*Tz;

		// Compute tangent
		Psi = gamma + beta*signum(TdStrain*Tz);
		Phi = 1.0 - pow(fabs(Tz),n)*Psi;
		double b1  = Ttheta*Phi;
		Ttangent = alpha*ko + (1-alpha)*ko*b1;
	}

	
	// Modify stress and tangent stiffness upon a reversal (refer to Drosos et al., 2012)
	if (mkur == 1) {
		if ( ((Tstress < CmaxStress) && (Tstress > 0.75*CmaxStress) && (Tstress < Cstress)) ||
			 ((Tstress > CminStress) && (Tstress < 0.75*CminStress) && (Tstress > Cstress)) ) {
			Ttangent = ko; //modify tangent
			dStress = ko*TdStrain;
			opserr << "Drosos: Tstress=" << Tstress << " CmaxStress=" << CmaxStress << " CminStress=" << CminStress << "\tTdStrain=" << TdStrain << " dStress=" << dStress; //forced debug
			opserr << " TmodStress=" << CmodStress + dStress;
			opserr << endln;
		}
		else dStress = Tstress - Cstress;
	}
	else dStress = Tstress - Cstress;
	
	TmodStress = CmodStress + dStress;


	// Eventually update the attained max and min strain
	if (Tstrain > CmaxStrain) TmaxStrain = Tstrain; else TmaxStrain = CmaxStrain;
	if (Tstrain < CminStrain) TminStrain = Tstrain; else TminStrain = CminStrain;

	// Eventually update the attained max and min stress
	if (Tstress > CmaxStress) TmaxStress = Tstress; else TmaxStress = CmaxStress;
	if (Tstress < CminStress) TminStress = Tstress; else TminStress = CminStress;


	//Define the variable recorded for debug
	TempVar = Ttangent;


    return 0;
}



double BoucWenGG::getStress(void) { return TmodStress; }

double BoucWenGG::getInitialTangent(void) { return ko; } //NON MI TORNA DIMENSIONALMENTE!!!!

double BoucWenGG::getTangent(void) { return Ttangent; }

double BoucWenGG::getStrain(void) { return Tstrain; }


// Add methods used by the class recorder
double BoucWenGG::getReferenceDuctility(void) { return Tmur; }

double BoucWenGG::getDeteriorationParameter(void) { return Ttheta; }

double BoucWenGG::getTempVar(void) { return TempVar; }



int BoucWenGG::commitState(void)
{
    // Commit trial history variables
    Cstrain = Tstrain;
	Cz = Tz;
	CmaxStrain = TmaxStrain;
	CmaxRevStrain = TmaxRevStrain;
	CminStrain = TminStrain;
	CminRevStrain = TminRevStrain;
	CdStrain = TdStrain;
	Cmurp = Tmurp;
	Cmurn = Tmurn;
	Cstress = Tstress;
	CmaxStress = TmaxStress;
	CminStress = TminStress;
	CmodStress = TmodStress;
	CNrevp = TNrevp;
	CNrevn = TNrevn;

    return 0;
}



int BoucWenGG::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}



int BoucWenGG::revertToStart(void)
{
    Tstrain = 0.0;
	Cstrain = 0.0;
	Tz = 0.0;
	Cz = 0.0;
	TmaxStrain = 0.0;
	CmaxStrain = 0.0;
	TmaxRevStrain = 0.0;
	CmaxRevStrain = 0.0;
	TminStrain = 0.0;
	CminStrain = 0.0;
	TminRevStrain = 0.0;
	CminRevStrain = 0.0;
	TdStrain = 0.0;
	CdStrain = 0.0;
	Tmurp = 0.0;
	Cmurp = 0.0;
	Tmurn = 0.0;
	Cmurn = 0.0;
	TNrevp = 0;
	CNrevp = 0;
	TNrevn = 0;
	CNrevn = 0;

	Tstress = 0.0;
	Cstress = 0.0;
	TmaxStress = 0.0;
	CmaxStress = 0.0;
	TminStress = 0.0;
	CminStress = 0.0;
	TmodStress = 0.0;
	CmodStress = 0.0;
	Ttangent = ko; //NON MI TORNA DIMENSIONALMENTE!!!!

	Tmur = 0.0;
	Ttheta = 0.0;
	TempVar = 0.0;

    return 0;
}



UniaxialMaterial* BoucWenGG::getCopy(void)
{
    BoucWenGG *theCopy =
	new BoucWenGG(this->getTag(), alpha, ko, strainY, n, gamma, beta, s1, s2, mkur, tolerance, maxNumIter);
    	
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->TmaxStrain = TmaxStrain;
    theCopy->TmaxRevStrain = TmaxRevStrain;
	theCopy->CmaxRevStrain = CmaxRevStrain;
	theCopy->TminStrain = TminStrain;
	theCopy->CminStrain = CminStrain;
	theCopy->TminRevStrain = TminRevStrain;
	theCopy->CminRevStrain = CminRevStrain;
	theCopy->TdStrain = TdStrain;
	theCopy->CdStrain = CdStrain;
	theCopy->Tmurp = Tmurp;
	theCopy->Cmurp = Cmurp;
	theCopy->Tmurn = Tmurn;
	theCopy->Cmurn = Cmurn;
	theCopy->CNrevp = CNrevp;
	theCopy->CNrevn = CNrevn;

    theCopy->Tstress = Tstress;
	theCopy->Cstress = Cstress;
	theCopy->TmaxStress = TmaxStress;
	theCopy->CmaxStress = CmaxStress;
	theCopy->TminStress = TminStress;
	theCopy->CminStress = CminStress;
	theCopy->TmodStress = TmodStress;
	theCopy->CmodStress = CmodStress;
    theCopy->Ttangent = Ttangent;

    return theCopy;
}



// Methods used in parallel processing with OpenSees
int BoucWenGG::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int BoucWenGG::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}



void BoucWenGG::Print(OPS_Stream &s, int flag)
{
    s << "BoucWenGG, tag: " << this->getTag() << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
	s << "  strainY: " << strainY << endln;
    s << "  n: " << n << endln;
    s << "  gamma: " << gamma << endln;
    s << "  beta: " << beta << endln;
    s << "  s1: " << s1 << endln;
    s << "  s2: " << s2 << endln;
	s << "  mkur: " << mkur << endln;
}



// Override the recorder
Response* BoucWenGG::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  theOutput.tag("BoucWenGGOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());

  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }

  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }

  // stress and strain
  else if ((strcmp(argv[0],"stressStrain") == 0)    || 
	       (strcmp(argv[0],"stressANDstrain") == 0) ||
	       (strcmp(argv[0],"stressAndStrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }

  // reference ductility
  else if (strcmp(argv[0],"ReferenceDuctility") == 0) {
    theOutput.tag("ResponseType", "mur");
    theResponse =  new MaterialResponse(this, 5, this->getReferenceDuctility());
  }

  // deterioration parameter
  else if (strcmp(argv[0],"DeteriorationParameter") == 0) {
    theOutput.tag("ResponseType", "theta");
    theResponse =  new MaterialResponse(this, 6, this->getDeteriorationParameter());
  }

  // stress, strain, ductility, and deterioration
  else if (strcmp(argv[0],"SSDD") == 0) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
	theOutput.tag("ResponseType", "mur");
	theOutput.tag("ResponseType", "theta");
    theResponse =  new MaterialResponse(this, 7, Vector(4));
  }

  // temp var
  else if (strcmp(argv[0],"TempVar") == 0) {
	theOutput.tag("ResponseType", "tempvar");
    theResponse =  new MaterialResponse(this, 8, this->getTempVar());
  }
	    
  theOutput.endTag();
  return theResponse;

}
 
int BoucWenGG::getResponse(int responseID, Information &matInfo)
{
  
	static Vector stressStrain(2);
	static Vector SSDD(4);

  switch (responseID) {
    case 1:
      matInfo.setDouble(this->getStress());
      return 0;
      
    case 2:
      matInfo.setDouble(this->getTangent());
      return 0;      

    case 3:
      matInfo.setDouble(this->getStrain());
      return 0;   

	case 4:
      stressStrain(0) = this->getStress();
      stressStrain(1) = this->getStrain();
      matInfo.setVector(stressStrain);
      return 0;

	case 5:
      matInfo.setDouble(this->getReferenceDuctility());
      return 0;   

	case 6:
      matInfo.setDouble(this->getDeteriorationParameter());
      return 0;    

    case 7:
      SSDD(0) = this->getStress();
      SSDD(1) = this->getStrain();
	  SSDD(2) = this->getReferenceDuctility();
	  SSDD(3) = this->getDeteriorationParameter();
      matInfo.setVector(SSDD);
      return 0;

	 case 8:
	  matInfo.setDouble(this->getTempVar());
      return 0;
      
  default:      
    return -1;

  }

}






