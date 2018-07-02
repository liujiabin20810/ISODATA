//=== File Prolog ================================================================
//      This code was developed by NASA, Goddard Space Flight Center, Code 588
//--- Contents -------------------------------------------------------------------
//      File 'Point.cc' contains the implementaion of 'Point' class as described in
//      'Point.h' in C++.
//
//--- Description ----------------------------------------------------------------
//  This calss was originally implemented for ISOCLUS algorithm but one may use it for
//  other image processing purposes.
//-- Notes:-----------------------------------------------------------------------
//
//-- Development History:--------------------------------------------------------
//   Date             Author                Reference
//   Description
//   
//   March 2002       Nargess Memarsadeghi  NASA GSFC, Code 588
//   Initial implementation
//
//--- DISCLAIMER---------------------------------------------------------------
//
//	This software is provided "as is" without any warranty of any kind, either
//	express, implied, or statutory, including, but not limited to, any
//	warranty that the software will conform to specification, any implied
//	warranties of merchantability, fitness for a particular purpose, and
//	freedom from infringement, and any warranty that the documentation will
//	conform to the program, or any warranty that the software will be error
//	free.
//
//	In no event shall NASA be liable for any damages, including, but not
//	limited to direct, indirect, special or consequential damages, arising out
//	of, resulting from, or in any way connected with this software, whether or
//	not based upon warranty, contract, tort or otherwise, whether or not
//	injury was sustained by persons or property or otherwise, and whether or
//	not loss was sustained from or arose out of the results of, or use of,
//	their software or services provided hereunder.
//--- Warning-------------------------------------------------------------------
//    This software is property of the National Aeronautics and Space
//    Administration.  Unauthorized use or duplication of this software is
//    strictly prohibited.  Authorized users are subject to the following
//    restrictions:
//    *   Neither the author, their corporation, nor NASA is responsible for
//        any consequence of the use of this software.
//    *   The origin of this software must not be misrepresented either by
//        explicit claim or by omission.
//    *   Altered versions of this software must be plainly marked as such.
//    *   This notice may not be removed or altered.
//
//=== End File Prolog=============================================================

#include <iostream>
#include <iomanip>
#include <new>
#include <math.h>
#include "Point.h"

using namespace std;
/****************************************************************************/
/* Constructor: Note it is assumed that the float* is dynamically allocated */
/* prior to be passed to constructor					    */
/****************************************************************************/
Point::Point(int dim, double* p)
{
	dimension = dim;
	point = p;
}
/****************************************************************************/
Point* Point::AllocPoint(int dim, double* p)
{
	dimension = dim;
	point = new double[dim];

	for (int i = 0; i < dim; i++)
		point[i] = p[i];

	return this;
}
/***************************************************************************/
/* Copy Constructor							   */
/**************************************************************************/
Point::Point(const Point& to_copy)
{
	dimension = to_copy.dimension;

	point = new double[dimension];
	if (!point)
	{
		cerr << "Memory Allocation Failed in Point::Point" << endl;
		cerr << "Exitting the program..." << endl;
		exit(1);
	}
	for (int i = 0; i < dimension; i++)
		point[i] = to_copy.point[i];

}

/****************************************************************************/
/* Default constructor							    */
/****************************************************************************/
Point::Point()
{
	dimension = 0;
	point = NULL;

}
/****************************************************************************/
/* Desstructor                                                              */
/****************************************************************************/
Point::~Point()
{
	if (point != NULL)
		delete[] point;

}
/******************************************************************************/
void Point::setPoint(double* p)
{
	point = p;

}
/******************************************************************************/
double Point::Norm2DistanceSquared(Point* p) const
{

	double distance = 0;

	if (dimension != p->getDimension())
	{
		cerr << "Points are not of the same dimension, cannot calculate their norm-2 distance" << endl;
		return -1;

	}

	//calculate the Eucladian distance between this point and the passed point.

	for (int i = 0; i < dimension; i++)
	{
		distance += pow((point[i] - p->point[i]), 2);

	}


	return distance;
}
/****************************************************************************/
void Point::print(ostream* out) const
{
	*out << "( ";
	for (int i = 0; i < dimension; i++)
	{
		if (i != 0)
			*out << ",";
		*out << setw(7) << point[i];

	}
	*out << ")\n";


}
/*************************************************************************/
int Point::getDimension() const
{
	return dimension;

}
/************************************************************************/
double* Point::getPoint()  const
{

	return point;

}
/************************************************************************/
int Point::operator==(const Point& p) const
{

	for (int i = 0; i < dimension; i++)
	{
		if (point[i] != p.point[i])
			return 0;
	}
	return 1;
}
/*************************************************************************/
/* addition of two points						 */
/*************************************************************************/
Point& Point::operator+(const Point& to_add)  const
{
	Point* result;

	if (dimension != to_add.dimension)
	{
		cout << "can't add points of different dimension:" << endl;
		this->print(&cout);
		to_add.print(&cout);
		exit(1);
	}
	try{
		double* r = new double[dimension];

		for (int i = 0; i < dimension; i++)
		{
			r[i] = point[i] + to_add.point[i];

		}

		result = new Point(dimension, r);
	}//try
	catch (bad_alloc exception){
		cout << "Exception occured in Point::operator+: " << exception.what() << endl;
		cout << "Exiting program..." << endl;
		exit(1);
	}


	return *result;
}
/*************************************************************************/
/* Multiplying a point by a number.					 */
/*************************************************************************/
Point& Point::operator*(double to_multiply)  const
{
	Point* result;

	try
	{
		double* r = new double[dimension];

		for (int i = 0; i < dimension; i++)
		{
			r[i] = point[i] * to_multiply;
		}

		result = new Point(dimension, r);

	} // try
	catch (bad_alloc exception)
	{
		cout << "Exception occured in Point::operator*: " << exception.what() << endl;
		cout << "Exiting program..." << endl;
		exit(1);
	}

	return *result;
}

/*************************************************************************/
/* Deviding a point by a number                                          */
/*************************************************************************/
Point& Point::operator/(double to_devide)  const
{

	if (to_devide == 0)
	{
		cout << " cannot devide a point by zero" << endl;
		exit(0);
	}

	double* r = new double[dimension];

	for (int i = 0; i < dimension; i++)
	{
		r[i] = point[i] / to_devide;

	}

	Point* result = new Point(dimension, r);

	return *result;
}
/************************************************************************/
/* adds a number to all coordinates of a point.			        */
/************************************************************************/
Point& Point::operator+(double to_add)  const
{
	double* r = new double[dimension];

	for (int i = 0; i < dimension; i++)
	{
		r[i] = point[i] + to_add;

	}

	Point* result = new Point(dimension, r);

	return *result;
}
/************************************************************************/
/* subtracts a number to all coordinates of a point.                         */
/************************************************************************/
Point& Point::operator-(double to_subtract)  const
{
	double* r = new double[dimension];

	for (int i = 0; i < dimension; i++)
	{
		r[i] = point[i] - to_subtract;

	}

	Point* result = new Point(dimension, r);

	return *result;
}


/*************************************************************************/

Point& Point::operator=(const Point& to_assign)
{
	int i;
	if (this != &to_assign)
	{
		if (dimension != to_assign.dimension)
		{
			dimension = to_assign.dimension;
			delete[] point;
			point = new double[dimension];
		}

		for (i = 0; i < dimension; i++)
			point[i] = to_assign.point[i];

	}

	return *this;
}
/***********************************************************************/
float Point::getCoordinate(int c)
{
	return point[c];

}
/**********************************************************************/
void Point::setCoordinate(int c, double value)
{

	point[c] = value;

}
