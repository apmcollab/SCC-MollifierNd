/*
 * NonSymmetricPolyMollifier1d.h
 *
 *  Created on: Dec 31, 2023
 *      Author: anderson
 */

#include <vector>
#include <iostream>

#ifndef  NON_SYMMETRIC_POLY_MOLLIFIER_1D_
#define  NON_SYMMETRIC_POLY_MOLLIFIER_1D_


#include "DoubleVectorNd/SCC_DoubleVector1d.h"

#include "LapackInterface/SCC_LapackMatrix.h"
#include "LapackInterface/SCC_LapackMatrixRoutines.h"

#include "OrthoFunctionNd/SCC_PolyFun.h"
#include "OrthoFunctionNd/SCC_OrthoPoly.h"

#include "SymmetricPolyMollifier1d.h"
#include "SplineFitNd/GaussLegendreTable.h"

#define MOLLIFIER_DEFAULT_ORDER       4
#define BASE_MOLLIFIER_DEFAULT_ORDER  4
#define DEFAULT_DIFFERENTIABLITY      6
/*
 *  Created on: Dec 31, 2023
 *      Author: anderson
 *
*/
//
//      NonSymmetricPolyMollifier1d
//
// A class that determines a continuous function such that
// the moments about xPos match up to the specified order
// (2 through 8) those of a delta function located at xPos.
//
// The function is non-zero only within [xMin,xMax] and
// and, if the differentiability parameter is set to p
// is p times continuously differentiable when extended
// to zero at the boundary where it vanishes.
//
//
//
// When xPos < (xMin+xMax)/2 function vanishes at xMax.
// (A left edge mollifier)
//
//         -
//      -     -
//   -           -
// -                    -
// |                           -
// |-------|-----------|---------------------|
// |      xPos                               |
// |                                         |
// xMin           (xMin+xMax)/2             xMax
//
//
// When xPos > (xMin+xMax)/2 function vanishes at xMin
// (A right edge mollifier.)
//
//
//                                   -
//                              -        -
//                       -                  -
//                  -                        -
// |           -
// |-------------------|--------------|------|
// |                                 xPos    |
// |                                         |
// xMin           (xMin+xMax)/2             xMax
//
/*
#############################################################################
#
# Copyright 2023 - Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

class NonSymmetricPolyMollifier1d
{
public:



NonSymmetricPolyMollifier1d()
{
	initialize();
}

NonSymmetricPolyMollifier1d(const NonSymmetricPolyMollifier1d& M)
{
	initialize(M);
}

NonSymmetricPolyMollifier1d(double xPos, double xMin, double xMax, double strength, int order, int differentiability)
{
	initialize(xPos, xMin, xMax, strength,order, differentiability);
}

void initialize()
{
	xMin = 0.0;
	xMax = 0.0;
    xPos = 0.0;

    strength = 0.0;

    order             = 0;
    differentiability = 0;
    baseOrder         = 0;
    nGauss            = 0;

	xGauss.clear();
	wtGauss.clear();

    baseMollifier.initialize();
    legendreFun.clear();

    basisCoeff.clear();
    basisFactors.clear();

    polyFactor.initialize();
    DpolyFactor.initialize();
	verboseFlag = false;

	gaussLegendreTable.initialize();
}

void initialize(const NonSymmetricPolyMollifier1d& M)
{
	initialize();
	xMin = M.xMin;
	xMax = M.xMax;
    xPos = M.xPos;

    strength = M.strength;

    order             = M.order;
    differentiability = M.differentiability;
    baseOrder         = M.baseOrder;
    nGauss            = M.nGauss;

	xGauss            = M.xGauss;
	wtGauss           = M.wtGauss;

    baseMollifier.initialize(M.baseMollifier);
    legendreFun  = M.legendreFun;
    basisCoeff   = M.basisCoeff;
    basisFactors = M.basisFactors;

    polyFactor.initialize(M.polyFactor);
    verboseFlag = M.verboseFlag;

    gaussLegendreTable.initialize();
}


void initialize(double xPos, double xMin, double xMax, double strength)
{
	initialize(xPos, xMin, xMax, strength,MOLLIFIER_DEFAULT_ORDER, DEFAULT_DIFFERENTIABLITY);
}


void setVerbose(bool val = true)
{
	verboseFlag = val;
}

void clearVerbose()
{
	verboseFlag = false;
}

void initialize(double xPos, double xMin, double xMax, double strength, int order, int differentiability)
{
    this->xMin     = xMin;
    this->xMax     = xMax;
    this->xPos     = xPos;
    this->strength = strength;

    this->order         = order;
    if(this->order > 8) {this->order = 8;}
	if(this->order < 1) {this->order = 1;}

    this->differentiability = differentiability;


    /////////////////////////////////////////////////////////////////////
    // Initialization of default base order
    //
    // Tricky part here when using a base mollifier of order 4 in order
    // to obtain reasonable mollifiers for high even order mollifiers and
    // low differentiability. Have to work out exactly why; now parameters
    // are being determined by experimentation.
    //
    // Consistent results for differentiability >= 4 and orders >= 4.
    //
    /////////////////////////////////////////////////////////////////////

    if(this->baseOrder == 0)
    {
    this->baseOrder  = BASE_MOLLIFIER_DEFAULT_ORDER;
    if((this->order < baseOrder)||((differentiability < 4) and (this->order%2 == 0)))
    {baseOrder = 2;}
    }

    if(this->order < baseOrder) {baseOrder = 2;}

    nGauss = 32;
	xGauss.resize(nGauss);
	wtGauss.resize(nGauss);

	double xMinStar = -1.0;
	double xMaxStar =  1.0;
	gaussLegendreTable.getLegendreNodesAndWeights(xMinStar,xMaxStar,nGauss,xGauss, wtGauss);

	createNormalizedMollifier();
}

void resetOrderAndDifferentiability(int order, int differentiability)
{
	this->order = order;
	if(this->order > 8) {this->order = 8;}
	if(this->order < 1) {this->order = 1;}

    this->differentiability = differentiability;

	if(this->differentiability > 9) {this->differentiability = 9;}
	if(this->differentiability < 1) {this->differentiability = 1;}

	createNormalizedMollifier();
}

void setBaseMollifierOrder(int val = BASE_MOLLIFIER_DEFAULT_ORDER)
{
	bool resetMollifierFlag = false;
	if(baseOrder != 0)
	{
		 resetMollifierFlag = true;
	}

    baseOrder = val;

    if(resetMollifierFlag)
    {
    createNormalizedMollifier();
    }
}

//
// The normalized mollifier evaluation is obtained by evaluating
//
// baseMollifer(x)*polyFactor(x)
//
// where polyFactor is the polynomial factor defined by the relation
//
//                                                         k = order-1
// baseMollifer(x)*polyFactor(x)  = baseMollifer(x)* (1 + SUM  basisFactor[k]*basisCoeff[k]*legendreFun[k](x));
//                                                         k = 0
//
// specifically
//                        k = order-1
// polyFactor(x)  = 1 + SUM  basisFactor[k]*basisCoeff[k]*legendreFun[k](x));
//                        k = 0
//
//

double operator()(double x) const
{
	if(order == 0) {return 0.0;}
	double xBar = -1.0 + 2.0*((x-xMin)/(xMax-xMin));

	if(xBar < -1.0) {return 0.0;}
	if(xBar >  1.0) {return 0.0;}

	return (baseMollifier(xBar)*polyFactor(xBar))/((xMax-xMin)*0.5);

	//
    // Alternate evaluation summing contributions individually
	/*
	double baseMollifierXbar = baseMollifier(xBar);
	double val   = basisFactors[0]*baseMollifierXbar;

	for(long k = 0; k < order; k++)
	{
		val += basisCoeff[k]*basisFactors[k]*baseMollifierXbar*legendreFun[k](xBar);
	}
    return val/((xMax-xMin)*0.5);
    */
}

	//  Returns a std::function that is bound to the evaluation operator of *this

std::function<double(double)> getEvaluationPtr() const
{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return F;
}


double derivative(double x) const
{
	if(order == 0) {return 0.0;}
	double xBar = -1.0 + 2.0*((x-xMin)/(xMax-xMin));
	if(xBar >  1.0) {return 0.0;}

	double val =
	baseMollifier.derivative(xBar)*polyFactor(xBar) + baseMollifier(xBar)*DpolyFactor(xBar);
	val *= 4.0/((xMax-xMin)*(xMax-xMin));

    return val;
}

void createNormalizedMollifier()
{

	double scaleFactor;

    double z =  -1.0 + 2.0*((xPos-xMin)/(xMax-xMin)); // Z = -1 => mollifier located at left edge
                                                      // Z =  1 => mollifier located at right edge;
    // Base the non-symmetric mollifier on a
    // shifted and scaled mollifier over [-1,1]

    double xCent       = z;
    double baseRadius  = std::max(1.0-z,1.0+z);

    baseMollifier.initialize(xCent, baseRadius, 1.0);

    // For order 0, just scale base mollifier to have unit integral

    if(order == 1)
    {
    	polyFactor.initialize(0);
    	baseMollifier.setOrder(2);
        baseMollifier.setDifferentiability(differentiability);
        scaleFactor = 0.0;
        for(long i = 0; i < nGauss; i++)
    	{
    		scaleFactor += baseMollifier(xGauss[i])*wtGauss[i];
    	}
        polyFactor[0] = 1.0/scaleFactor;
        DpolyFactor.initialize(0);
        return;
    }

    baseMollifier.setOrder(baseOrder);
    baseMollifier.setDifferentiability(differentiability);

    // Create Legendre functions over [-1,1] for the
    // polynomial function multipliers of the baseMollifier
    // that are used as a basis for the mollifer

    SCC::OrthoPoly orthogonalPoly(SCC::OrthoPoly::Legendre);
    legendreFun = orthogonalPoly.getOrthoPolyArray(order-1);

    // Create product basis functions values

    std::vector<SCC::DoubleVector1d> basisFun(order);
    basisFactors.resize(order);


    for(long k = 0; k < order; k++)
    {
    	basisFun[k].initialize(nGauss);
    	for(long i = 0; i < nGauss; i++)
    	{
    		basisFun[k](i) = baseMollifier(xGauss[i])*legendreFun[k](xGauss[i]);
    	}

    	// L2 normalize

    	scaleFactor = 0.0;
    	for(long i = 0; i < nGauss; i++)
    	{
    		scaleFactor += basisFun[k](i)*basisFun[k](i)*wtGauss[i];
    	}
    	scaleFactor = std::sqrt(1.0/scaleFactor);
    	basisFactors[k] = scaleFactor;
    	for(long i = 0; i < nGauss; i++)
    	{
    		basisFun[k](i) *= scaleFactor;
    	}
    }

    // Construct moment equations

    double integralVal;

    int mx;
    std::function<double(double)> momentFunction = [&mx,&xCent](double x)
	{
    	double val = std::pow((x-xCent),mx);
    	return val;
	};

    SCC::DoubleVector1d momentValues(nGauss);

    SCC::LapackMatrix Rstar(order,1);
    double integralStar;

    if(verboseFlag)
    {
    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Base mollifier correction moments " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    }

	for(mx = 0; mx < order; mx++)
	{
		for(long i = 0; i < nGauss; i++)
		{momentValues(i) = momentFunction(xGauss[i]);}

		integralStar = 0.0;
		for(long i = 0; i < nGauss; i++)
		{
			integralStar += momentValues(i)*basisFactors[0]*baseMollifier(xGauss[i])*wtGauss[i];
		}

		if(mx == 0)
		{
			Rstar(mx) = -(integralStar - 1.0) ;
		}
		else
		{
			Rstar(mx) = -integralStar;
		}
		if(verboseFlag) {std::cout << " mx : " << Rstar(mx) << std::endl;}
	}

	SCC::LapackMatrix Mstar(order,order);

	for(int j = 0; j < order; j++)
	{
		for(mx = 0; mx < order; mx++)
		{
			for(long p = 0; p < nGauss; p++)
			{momentValues(p) = momentFunction(xGauss[p]);}

			integralVal = 0.0;
			for(long p = 0; p < nGauss; p++)
			{
				integralVal += momentValues(p)*basisFun[j](p)*wtGauss[p];
			}

			Mstar(mx,j) = integralVal;
		}
    }

    // M.printDense(std::cout, precision);

    SCC::DGESVD dgesvd;
    dgesvd.computeSVD(Mstar);

    if(verboseFlag)
    {
    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Correction matrix singular values " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    for(long i = 0; i < order; i++)
    {
    std::cout << dgesvd.singularValues[i] << std::endl;
    }

    std::cout << "\nCondition number " << dgesvd.singularValues[0]/dgesvd.singularValues[order-1] << std::endl;
    }

    // Set up and solve for the coefficients

    SCC::LapackMatrix Cstar(order,order);

    for(long i = 0; i < order; i++)
    {
    for(long j = 0; j < order; j++)
    {
    Cstar(i,j) = Mstar(i,j);
    }}

    //Cstar.printDense(std::cout, 5);

    SCC::DGESVX               dgesvx;
    dgesvx.applyInverse(Cstar, Rstar);

    // Capture the coefficients for the normalized mollifier

    basisCoeff.resize(order);
    for(long i = 0; i < order; i++)
    {
    	basisCoeff[i] = Rstar(i);
    }

    polyFactor.initialize(order-1);
    polyFactor[0] = basisFactors[0];
	for(long k = 0; k < order; k++)
	{
		polyFactor = polyFactor + basisCoeff[k]*basisFactors[k]*legendreFun[k];
	}

	DpolyFactor = polyFactor.differentiate();

    if(verboseFlag)
    {
    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Normalized mollifier polynomial factor  " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout.precision(17);
    std::cout << polyFactor << std::endl;

    /*
    for(long i = 0; i < order; i++)
    {
    	if(i == 0)
    	{
    	std::cout << i << " : Coeff = " << (1 + basisCoeff[i])*basisFactors[i] << std::endl;
    	}
    	else
    	{
    	std::cout << i << " : Coeff = " << basisCoeff[i]*basisFactors[i] << std::endl;
    	}
    }
    */
    }

    }

    void momentCheck()
    {
    double xCent = xPos;
    int    mx;

    std::function<double(double)> momentFunction = [&mx,&xCent]
    (double x)
    {
        double val = std::pow((x-xCent),mx);
    	return val;
    };

    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Non-Symmetric Mollifier moment error " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    double integralStar;
    double momentError = 0.0;
    for(mx = 0; mx < order; mx++)
    {
    	integralStar = 0.0;
    	for(long i = 0; i < nGauss; i++)
    	{integralStar += momentFunction(xGauss[i])*this->operator()(xGauss[i])*wtGauss[i];}

    	if(mx == 0)
    	{
    		momentError = std::max(momentError,std::abs(integralStar - 1.0));
    		std::cout << " mx_" << mx << "_Error : " << std::abs(integralStar - 1.0) << std::endl;
        }
    	else
    	{
    		momentError = std::max(momentError,std::abs(integralStar));
    		std::cout << " mx_" << mx << "_Error : " <<  std::abs(integralStar) << std::endl;
    	}
    }

    }



    double  xMin;    // The mollifier is non-zero over the interval [xMin,xMax]
	double  xMax;
    double  xPos;    // The location mollifier within the interval [xMin,xMax]

    double    strength;    // The strength of the mollifier

    int             order;    // Order of the mollifier = 2-8 (odd order specification results in next higher order)
    int differentiability;    // The differentiability of the mollifier (1-9)
    int         baseOrder;    // The order of the base mollifier
    int            nGauss;    // The number of Gaussian quadrature nodes used

    GaussLegendreTable gaussLegendreTable;
    std::vector<double>            xGauss;
    std::vector<double>           wtGauss;

    SymmetricPolyMollifier1d baseMollifier;
    std::vector<SCC::PolyFun>  legendreFun;

    std::vector<double> basisCoeff;
    std::vector<double> basisFactors; // L2 normalization factors

    SCC::PolyFun        polyFactor;   // The polynomial factor multiplying the base mollifier
    SCC::PolyFun       DpolyFactor;   // The derivative of the polynomial factor
    bool verboseFlag;
};

#undef MOLLIFIER_DEFAULT_ORDER
#undef BASE_MOLLIFIER_DEFAULT_ORDER
#undef DEFAULT_DIFFERENTIABLITY

#endif //NON_SYMMETRIC_POLY_MOLLIFIER_1D_



