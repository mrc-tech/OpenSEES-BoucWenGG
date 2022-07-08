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
                                                                        
// $Revision: 0.0
// $Date: 2013/10/11
// $Source: 
                                                                        
#ifndef BoucWenGG_h
#define BoucWenGG_h

// Written:  Andrea Lucchini (andrea.lucchini@uniroma1.it)
// Modified: Andrea Marchi   (andrea.marchi@uniroma1.it) [19/02/2021]
//
// Description: This file contains the class definition for the
// BoucWenGG Uniaxial Material (refer to Gerolymos and Gazetas, 2005, and Drosos et al., 2012).
// BoucWenGG is a modified version of the classic Bouc-Wen model that includes also cyclic deterioration. 
//
// What: 

#include <UniaxialMaterial.h>

class BoucWenGG : public UniaxialMaterial
{
  public:
    BoucWenGG(int tag, 
		    double alpha,
		    double ko,
			double strainY,
		    double n,
		    double gamma,
		    double beta,
		    double s1,
		    double s2,
			int mkur, 
		    double tolerance,
		    int maxNumIter);  
    BoucWenGG();    

    ~BoucWenGG();

	const char *getClassType(void) const {return "BoucWenGG";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

	// Methods called by the class recorder
	double getReferenceDuctility(void);		
	double getDeteriorationParameter(void);
	double getTempVar(void); // for debug 

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);  

	double signum(double); //sign function [sgn(x)]

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
	
	// New recorder for this class
	Response *setResponse (const char **argv, int argc, OPS_Stream &theOutputStream);
    int getResponse (int responseID, Information &matInformation);  

  protected:
    
  private:

  // Material parameters
    double alpha;	// final stiffness / initial stifness
    double ko;		// initial stiffness
	double strainY;	// yielding strain (or deformation)
	double n;
    double gamma;
    double beta;
    double s1;
    double s2;
	int mkur;		// if mkur is == 1 -> the stiffness value modification upon reversal (Drosos et al., 2012) is applied

    // History variables (trial and commited)
    double Tstrain, Cstrain;
	double Tz, Cz;

	double TmaxStrain, CmaxStrain;			// attained max strain (positive value)
	double TmaxRevStrain, CmaxRevStrain;	// attained max strain at reversal
	double TminStrain, CminStrain;			// attained min strain (negative value)
	double TminRevStrain, CminRevStrain;	// attained min strain at reversal
	double TdStrain, CdStrain;				// strain increment 
	double Tmurp, Cmurp;					// reference strain ductility (Drosos et al., 2012), positive value
	double Tmurn, Cmurn;					// reference strain ductility (Drosos et al., 2012), negative value
	double Tstress, Cstress;				// stress
	double TmaxStress, CmaxStress;			// attained max stress (positive value)
	double TminStress, CminStress;			// attained min stress (positive value)
	double TmodStress, CmodStress;			// modified stress upon a reversal (Drosos et al., 2012)
	int TNrevp, CNrevp;						// counter of reversals in the positive stress domain
	int TNrevn, CNrevn;						// counter of reversals in the negative stress domain

	// Ohter variables
	double Ttangent;						// tangent stiffness
	double Tmur;							// reference strain ductility
	double Ttheta;							// deterioration parameter
	double TempVar;							// variable recorded for debug 
	double tolerance;						// tolerance in the calculation of the internal variable z of the model (an iterative Newton scheme is applied)
	int maxNumIter;							// maximum number of iterations in the calculation of the internal variable z


	// Auxiliary functions
	double fun(double z, double z_old); //function to be evaluated internally (fun(z,z_old) == 0) do be solved for z
	

};


#endif



