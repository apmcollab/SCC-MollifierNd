#include <vector>

#ifndef MOMENT_INDICES_
#define MOMENT_INDICES_
//
// A class with a single utility function for determining all 
// moment indices of a specified order. 
//

/*
    // Code snippet demonstration of the listing af all 3d moments

    long xIndex; long yIndex; long zIndex;
	MomentIndices momentIndices;

    int dimension = 3;
    std::vector<std::vector < int > > indices;

    momentIndices.getMomentIndices(order,dimension,indices);

	for(long k = 0; k < indices.size(); k++)
	{
    xIndex = indices[k][0];
    yIndex = indices[k][1];
    zIndex = indices[k][2];
    cout << xIndex << " " << yIndex << " " << zIndex << endl;
	}

 */

 /*
#############################################################################
#
# Copyright 2015-16-16 Chris Anderson
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

class MomentIndices
{

	public:

	void getMomentIndices(int order, int dimension, std::vector<std::vector <int> >& indices)
    {
    std::vector<int> indexValues(dimension,0);
    int coordinate = 0;

	getMomentIndex(order, coordinate, indexValues,indices);
    }

//
// Recursive method to find all indices of moments up to a given order.
// The dimension of the coordinate system is determined by the
// size of the indexValue std::vector.
//
//  Initial inputs :
//
//  order      = order of the moments
//  coordinate = 0
//  indexValues = std::vector of int's of length equal to the dimension of the coordinate system
//
//
// Output :

// indices : an array of indices = an std::vector of vectors containing the indices associated
//           with the specified order

void getMomentIndex(int order, int coordinate, std::vector<int> indexValues,
std::vector<std::vector <int> >& indices)
{
    if(order == 0)
    {
    	indices.push_back(indexValues);
    	return;
    }
    if(coordinate == (int)indexValues.size() - 1)
    {
    	indexValues[coordinate] = order;
    	indices.push_back(indexValues);
    	return;
    }

    int sub_order;
    for(int i = 0;  i <= order; i++)
    {
    indexValues[coordinate] = i;
    sub_order  = order - i;
    if(coordinate < (int)indexValues.size()-1) {getMomentIndex(sub_order, coordinate+1,indexValues,indices);}
    }
}
};

#endif
