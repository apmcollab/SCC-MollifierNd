/*
 * DiscretePolyMollifier1d.h
 *
 *  Created on: Jan 4, 2024
 *      Author: anderson
 */
 /*
 * A class that determines discrete grid values such that
 * the moments about xPos computed with the trapezoidal
 * method integration formula match up to the specified order
 * (2 through 8) those of a delta function located at xPos.
 *
 * The discrete mollifiers are obtained by sampling either 
 * the symmetric polynomial mollifier or non-symmetric 
 * polymollifier and then correcting the values to obtain
 * the desired moments. 
 *
 * In order for the correction equations to be well determined
 * the panel width - the number of panels determining the 
 * grid points on either side of xPos where values are distributed
 * must be greater than or equal to (order/2) + 1.
 *
 * If the grid points over which the values are distributed extends
 * past the endpoint of the interval, then the distribution is based
 * upon the non-symmetric polynomial mollifier, otherwise the 
 * distribution is based upon the symmetric mollifier.  
 *
 *
 */
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

#include <iostream>
#include <stdexcept>

#ifndef DISCRETE_POLY_MOLLIFIER_1D_
#define DISCRETE_POLY_MOLLIFIER_1D_

#include "GridFunctionNd/SCC_GridFunction1d.h"

#include "SymmetricPolyMollifier1d.h"
#include "NonSymmetricPolyMollifier1d.h"

#include "LapackInterface/SCC_LapackMatrix.h"
#include "LapackInterface/SCC_LapackMatrixRoutines.h"


#define MOLLIFIER_DEFAULT_ORDER  4
#define DEFAULT_DIFFERENTIABLITY 6


class DiscretePolyMollifier1d
{
public:

	DiscretePolyMollifier1d()
	{
		initialize();
	}

	DiscretePolyMollifier1d(const DiscretePolyMollifier1d& D)
	{
	    order             = D.order;
		differentiability = D.differentiability;
		verboseFlag       = D.verboseFlag;
		periodicFlag      = D.periodicFlag;
		symmetricFlag     = D.symmetricFlag;
	}

	void initialize()
	{
	    order             = MOLLIFIER_DEFAULT_ORDER;
		differentiability = DEFAULT_DIFFERENTIABLITY;
		verboseFlag       = false;
		periodicFlag      = false;
		symmetricFlag     = true;
		symmetricPolyMollifier.initialize();
		nonSymmetricPolyMollifier.initialize();
	}

	void initialize(const DiscretePolyMollifier1d& D)
	{
	    order             = D.order;
		differentiability = D.differentiability;
		verboseFlag       = D.verboseFlag;
		periodicFlag      = D.periodicFlag;
		symmetricFlag     = D.symmetricFlag;
		symmetricPolyMollifier.initialize(D.symmetricPolyMollifier);
		nonSymmetricPolyMollifier.initialize(D.nonSymmetricPolyMollifier);
	}


    void setOrder(int order)
	{
	this->order = order;
	}

    void setDifferentiability(int diffOrder)
	{
	differentiability = diffOrder;
	}

	int getDifferentiablity()
	{
	return differentiability;
	}

	int getOrder() const
	{return this->order;}

    void setVerbose(bool val = true)
	{
		verboseFlag = val;
	}

	void clearVerbose()
	{
		verboseFlag = false;
	}

    void createDiscreteMollifier(double xPos, long panelWidth, double strength, SCC::GridFunction1d& B)
	{
    	double xMin  = B.getXmin();
    	double xMax  = B.getXmax();
    	long xPanels = B.getXpanelCount();
    	double hx    = B.getHx();

        long xGridIndexBase = std::round((xPos - xMin)/hx);
	    long xGridIndxMin   = (xGridIndexBase - panelWidth);
		long xGridIndxMax   = (xGridIndexBase + panelWidth);

		double tol     = 1.0e-14;
		double tolStar = tol*(1.0 + std::abs(xPos));


		if(xPos > (xGridIndexBase*hx+xMin) + tolStar) {xGridIndxMax += 1;}
	    if(xPos < (xGridIndexBase*hx+xMin) - tolStar) {xGridIndxMin -= 1;}

	    if (xPanels < 2*panelWidth)
		{
			std::string errMsg = "\n  DiscretePolyMollifier: Insufficient panels in GridFunction1d. \n";
			errMsg            += "                         GridFunction1d panels must be >= 2*panelWidth \n";
			errMsg            += "                         Input panels         = " + std::to_string(xPanels) + "\n";
			errMsg            += "                         Input panel width    = " + std::to_string(panelWidth) + "\n";
			errMsg            += "                         Minimum panels       = " + std::to_string(2*panelWidth) + "\n";
			throw std::runtime_error(errMsg);
		}

	    if (panelWidth < (order/2) + 1)
		{
			std::string errMsg = "\n  DiscretePolyMollifier: Insufficient discrete mollifier panel width. \n";
			errMsg            += "                         Mollifier panel width must be >= (order/2) + 1 \n";
			errMsg            += "                         Input order          = " + std::to_string(order) + "\n";
			errMsg            += "                         Input panel width    = " + std::to_string(panelWidth) + "\n";
			errMsg            += "                         Minimum panel width  = " + std::to_string(order/2  + 1) + "\n";
			throw std::runtime_error(errMsg);
		}

    	if((((xPos - panelWidth*hx) >= xMin)&&((xPos + panelWidth*hx) <= xMax))||(periodicFlag))
    	{
    		symmetricFlag = true;
    		if(verboseFlag)
    		{std::cout << "Symmetric mollifier construction" << std::endl;}
    	}
    	else
    	{
    		symmetricFlag = false;
    		if(verboseFlag)
    		{std::cout << "Non-symmetric mollifier construction" << std::endl;}
    	}

    	// Create local mollifier scaled to [-1,1]

        double  localXmin    = -1.0;
        double  localXmax    =  1.0;

        long    localXpanels;
        double  localHx;
        double  hxFactor;
        double  localXpos;

        if(symmetricFlag)
        {
        	localXpanels   = xGridIndxMax - xGridIndxMin;
        	localXpos      = -1.0 + 2.0*(xPos - (xGridIndxMin*hx+xMin))/(localXpanels*hx);
        	localHx        = (localXmax-localXmin)/localXpanels;
        	hxFactor       = (hx/localHx);
        	symmetricPolyMollifier.initialize(localXpos, panelWidth*localHx,strength);
        	symmetricPolyMollifier.setOrder(order);
        	symmetricPolyMollifier.setDifferentiability(differentiability);
        }
        else
        {
        	localXpanels = 2*panelWidth;
        	localHx      = (localXmax-localXmin)/localXpanels;
        	hxFactor     = (hx/localHx);

        	if((xPos - panelWidth*hx) < xMin)  // mollifier located near left edge
        	{
        		localXpos =  -1.0 + 2.0*(xPos - xMin)/(localXpanels*hx);
        	}
        	else                              // mollifier located near right edge
        	{
        		localXpos =   1.0 - 2.0*(xMax - xPos)/(localXpanels*hx);
        	}
        	nonSymmetricPolyMollifier.initialize(localXpos, localXmin, localXmax, strength, order, differentiability);
        }

        if(verboseFlag)
        {
        std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        std::cout << "Discrete source construction " << std::endl;
        std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        }

        double integralVal;

        int mx;

        std::function<double(double)> momentFunction = [&mx,&localXpos](double x)
    	{
        double val = std::pow((x-localXpos),mx);
    	return val;
    	};

        SCC::GridFunction1d FA(localXpanels,localXmin,localXmax);
        SCC::GridFunction1d xGridFun(localXpanels,localXmin,localXmax);

        if(symmetricFlag){xGridFun.specify(symmetricPolyMollifier);}
        else             {xGridFun.specify(nonSymmetricPolyMollifier);}

        double scaleFactor = xGridFun.integralTrapezoidal()*hxFactor;
        if( scaleFactor > 0.1){xGridFun  *= 1.0/scaleFactor;}

        std::ios_base::fmtflags ff;
        int precision;
        int precisionCache;

        // Compute moments

        if(verboseFlag)
        {
        std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        std::cout << "Base mollifier discrete moments " << std::endl;
        std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;


        precision  = 3;
        ff         =  std::cout.flags();
        precisionCache = std::cout.precision(precision);
        for(mx = 0; mx < order+1; mx++)
        {
    	FA.specify(momentFunction);
    	FA *= xGridFun;
    	integralVal = FA.integralTrapezoidal()*hxFactor;
    	if(mx == 0)
    	{std::cout << " mx_" << mx << "_Error : " <<  std::scientific <<  std::right << std::setw(precision+7) << std::abs(integralVal  - 1.0) << " ";}
    	else
    	{std::cout << " mx_" << mx << "_Error : "  <<  std::right << std::setw(precision+7) << std::abs(integralVal) << " ";}
        std::cout << std::endl;
        }

        std::cout.flags(ff);
        std::cout.precision(precisionCache);
        }

        // Create Legendre basis [-1,1]

        SCC::OrthoPoly orthogonalPoly(SCC::OrthoPoly::Legendre);
        std::vector<SCC::PolyFun> legendreFun = orthogonalPoly.getOrthoPolyArray(order-1);

        SCC::GridFunction1d FB(localXpanels,localXmin,localXmax);

        precision = 5;
        ff =  std::cout.flags();
        precisionCache = std::cout.precision(precision);

        SCC::LapackMatrix M(order,order);

        for(int i = 0; i < order; i++)
        {
    	FA.specify(legendreFun[i].getEvaluationPtr());
    	FA *= xGridFun;

    	scaleFactor = FA.norm2()*std::sqrt(hxFactor);
    	FA *= 1.0/std::sqrt(scaleFactor);

    	for(mx = 0; mx < order; mx++)
    	{
    	FB.specify(momentFunction);
    	FB *= FA;

    	integralVal = FB.integralTrapezoidal()*hxFactor;
    	M(mx,i) = integralVal;
    	}
    }
    std::cout.flags(ff);
    std::cout.precision(precisionCache);

    dgesvd.computeSVD(M);

    if(verboseFlag)
    {
    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Discrete correction matrix singular values " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    for(long i = 0; i < order; i++)
	{
    	std::cout << dgesvd.singularValues[i] << std::endl;
	}

	std::cout << "\nCondition number " << dgesvd.singularValues[0]/dgesvd.singularValues[order-1] << std::endl;
    }

  //
    // Form right hand side
    //
    SCC::LapackMatrix R(order,1);
    for(mx = 0; mx < order; mx++)
    {
    	FA.specify(momentFunction);
    	FA *= xGridFun;
    	integralVal = FA.integralTrapezoidal()*hxFactor;

    	if(mx == 0)
    	{
    	R(mx) = -(integralVal - 1.0) ;
    	}
    	else
    	{
    	R(mx) = -integralVal;
    	}
    }

    // Form matrix

    SCC::LapackMatrix C(order,order);
    for(long i = 0; i < order; i++)
    {
    for(long j = 0; j < order; j++)
    {
    	C(i,j) = M(i,j);
    }}

    dgesvx.applyInverse(C, R);

    if(verboseFlag)
    {
    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Discrete correction coefficients " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    for(long i = 0; i < order; i++)
    {
    	std::cout << "Coeff " << R(i) << std::endl;
    }
    }
    // Add in correction

    SCC::GridFunction1d xGridFunStar(localXpanels,localXmin,localXmax);
    xGridFunStar = xGridFun;

    for(int i = 0; i < order; i++)
    {
    	FA.specify(legendreFun[i].getEvaluationPtr());

    	FA *= xGridFun;
    	scaleFactor = FA.norm2()*std::sqrt(hxFactor);
    	FA *= 1.0/std::sqrt(scaleFactor);
    	FA *= R(i);

    	xGridFunStar += FA;
    }

    if(verboseFlag)
    {
    SCC::GridFunction1d diffGrid(localXpanels,localXmin,localXmax);

    diffGrid  = xGridFunStar;
    diffGrid -= xGridFun;
    std::cout << std::endl;
    std::cout << "Correction_Norm_Rel_Inf " << diffGrid.normInf()/xGridFun.normInf() << std::endl;
    std::cout << "Correction_Norm_Rel_L2  " << diffGrid.norm2()/xGridFun.norm2() << std::endl;
    std::cout << std::endl;
    }

    B.setToValue(0.0);

    long periodicIndex;
	long Bpanels = B.getXpanelCount();

    if(symmetricFlag)
    {
		if(not periodicFlag)
		{
			for(long i = 1; i < localXpanels; i++)
			{
			B(xGridIndxMin + i) = xGridFunStar(i);
			}
		}
		else
		{
		    for(long i = 0; i <= localXpanels; i++)
			{
		    	if(xGridIndxMin + i >= Bpanels)
		    	{
		    	periodicIndex = (xGridIndxMin - Bpanels) + i;
		    	B(periodicIndex) = xGridFunStar(i);
		    	if(periodicIndex == 0){B(Bpanels) = xGridFunStar(i);}
		    	}
		    	else if (xGridIndxMin + i <= 0)
		    	{
		    	periodicIndex = (xGridIndxMin + Bpanels) + i;
		    	B(periodicIndex) = xGridFunStar(i);
		    	if(periodicIndex == Bpanels){B(0) = xGridFunStar(i);}
		    	}
		    	else
		    	{
		    	B(xGridIndxMin + i) = xGridFunStar(i);
		    	}
			}
		}
    }
    else
    {
    //
    // Propagate source values to grid, scaling appropriately
    //
    if((xPos - panelWidth*hx) < xMin)  // mollifier located near left edge
	{
    	for(long i = 0; i <= localXpanels; i++)
    	{
    		B(i) = xGridFunStar(i);
    	}
	}
    else         // mollifier located near right edge
	{
    	for(long i = 0; i <= localXpanels; i++)
    	{
    		B(xPanels-i) = xGridFunStar(localXpanels-i);
    	}
	}}

	}


	double periodicShift(double xMin, double xMax, double& x)
	{
		long shift;
		if(x < xMin)
		{
			shift = std::floor(std::abs(x-xMin)/(xMax-xMin)) + 1;
			x    += shift*(xMax-xMin);
		}
		else
		{
			shift = std::floor(std::abs(x-xMin)/(xMax-xMin));
			x    -= shift*(xMax-xMin);
		}
		return x;
	}

    void createPeriodicDiscreteMollifier(double xCent, long panelWidth, double strength, SCC::GridFunction1d& B)
	{
    	double xMin         = B.getXmin();
    	double xMax         = B.getXmax();
    	periodicFlag        = true;
    	periodicShift(xMin, xMax,  xCent);
    	createDiscreteMollifier(xCent, panelWidth, strength, B);
    	periodicFlag       = false;
	}

    void checkOrder(SCC::GridFunction1d& M, double xBar)
	{
	double xCent = xBar;
    int mx;
    std::function<double(double)> momentFunction = [&mx,&xCent](double x)
    {
        double val = std::pow((x-xCent),mx);
    	return val;
    };

    std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "Discrete mollifier moment error " << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;


    std::ios_base::fmtflags ff;
    int precision;
    int precisionCache;

    //
    // Compute moments
    //
    SCC::GridFunction1d FA(M);
    FA.setToValue(0.0);

    double integralVal;

    int polyMaxIndex = 10;

    precision  = 3;
    ff =  std::cout.flags();
    precisionCache = std::cout.precision(precision);
    for(mx = 0; mx < polyMaxIndex; mx++)
    {
    	FA.specify(momentFunction);
    	FA *= M;
    	integralVal = FA.integralTrapezoidal();
    	if(mx == 0)
    	{std::cout << " mx_" << mx << "_Error : " <<  std::scientific <<  std::right << std::setw(precision+7) << std::abs(integralVal  - 1.0) << " ";}
    	else
    	{std::cout << " mx_" << mx << "_Error : " << std::scientific <<  std::right << std::setw(precision+7) << std::abs(integralVal) << " ";}
        std::cout << std::endl;
    }

    std::cout.flags(ff);
    std::cout.precision(precisionCache);
	}


	long order;
	long differentiability;
    bool verboseFlag;
    bool periodicFlag;
    bool symmetricFlag;

    SymmetricPolyMollifier1d     symmetricPolyMollifier;
	NonSymmetricPolyMollifier1d  nonSymmetricPolyMollifier;

	SCC::DGESVD  dgesvd;
	SCC::DGESVX  dgesvx;

};

#undef MOLLIFIER_DEFAULT_ORDER
#undef DEFAULT_DIFFERENTIABLITY

#endif // DISCRETE_POLY_MOLLIFIER_1D_
