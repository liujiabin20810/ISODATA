//=== File Prolog ================================================================
//      This code was developed by NASA, Goddard Space Flight Center, Code 588
//      and computer science department of university of Maryland at College Park.
//--- Contents -------------------------------------------------------------------
//      File 'Image.h' contains the definition of 'Image' class in C++.
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

#ifndef IMAGE_H
#define IMAGE_H

#include "KM_ANN.h"
#include "Point.h"
#include <string>
#include <vector>

extern bool is_rand;
extern int  sample_seed;
/*******************************************************************************************/
/* PairDistanceNode contains the pair wise distance between two cluster centers c1, and c2.*/
/* we will need to save the distance as well as the clusters number of which c1, and c2 are*/
/* their centers.									   */
/*******************************************************************************************/
struct PairDistanceNode {

double dist;
int c1;
int c2;
 
 bool operator<(PairDistanceNode to_compare)  const                  
 {
 return  (dist < to_compare.dist);

 }

 bool operator<(double to_compare) const
 {
 return (dist < to_compare);
    
 } 

};

/**********************************************************************************************/
/* class 'image' supports structures for storing and displaying an image as well as performing*/
/*  the required image processing functionalities and structures for ISOCLUS algorithm.       */
/* ie: the main goal of this class is supporting the implementation of ISOCLUS algorithm.     */
/* one could use this class to implement other clustering algorithms.                         */
/**********************************************************************************************/  
class Image{

public:
	Image();
	Image(int row, int col, int bands, int NumClus, int SAMPRM);
	~Image();

	void SetNumClusters();
	void readImage(std::string* );
	void writeClassifiedImage(std::string);
	void CalculateDistances();
	void UpdateCenters();

        /* one can obtain a point of the image either by passing a number corresponding to */
        /* the pixel's location in the image from the begining (pixel_row*NumRow+pixel_col),*/
        /* or by passing pixel's coordinates (pixel_row, pixel_col)                         */
        /* Note when refering to a pixel by its coordinates, the first pixle of the image is*/
        /* refered to as (0,0)                                                              */
    
        Point* getPoint(int PointCount);
        Point* getPoint(int row, int col);
	void   setPoints(KMpointArray all);
	int    size();
	bool   WasDeleted();

	//Based on the desired NumClusters randomly select a set of sample centers.
         void   sampleCenters(); // parameter is true if data was generated , 
	                            // and false if it was read from file.

	//Select a set of initial sample indices of points from allPoints array.
	void samplePoints(double);

	// add a small purturbation
        void addNoise();

	//Distribute sampled data points among the present cluster centers.
	void   PutInCluster();

	//If any cluster's size is less than NumSamples delete that center and cluster.
	void   PostAnalyzeClusters();

	//Calculate the average distance of points in a cluster from their centers.
	void   CalculateAverageDistances();

	//Calculate the overall average distance of points from their cluster center.
	double  OverallAverageDistances();

	//Returns average of sum of squared distances of points from their cluster centers.
	double getDistortions();
 
	//Calculate the standard deviation vector (a point for each cluster).
	void CalculateSTDVector();

	//Find the maximum element of each column of STDVector.
	void CalculateVmax();

	//Decide if we need to split any particular clusters or not. Return a vector of
        //integers which contains the cluster numbers that need to split.
        std::vector<int> ShouldSplit(double stdv);


	//splits clusters whose index appear in the passed vector.
	void Split(std::vector<int> to_split);

	//computes the pairwise distance between cluster centers
	void ComputeCenterDistances(); 

	std::vector<PairDistanceNode> FindLumpCandidates(double lump, int MAXPAIR);
        void Lump(const std::vector<PairDistanceNode>& to_lump);     

	void printCoordinates(int pos);
	int  getNumCenters();
        void preFinalClustering();
        void generateReport();

		void getClusterLabel(int* label);
private:
	// Each Image is a two dimensional array of unsigned chars (8 bits per pixel).
	// first dimension is of size NumRows*NumColumns of the image (all pixels) and
	// second dimension is of size NumBands. Thus, row one of the 'image' array 
	// will have all pixels in band one, row two will have all pixels in band 2, etc. 

	unsigned char** io_image;

	//
	int NumBands;
	int NumRow;
	int NumCol;
	int ImageSizeInByte; //which is basically NumRow*NumCol
	int NumClusters;    //desired number of clusters to get after classification
        unsigned int NumSamples;
	bool Deleted;

	Point** allPoints;
	std::vector<Point*> centers;
	//indices of sampled points selected from allPoints array.
	std::vector<int> samples;

	//we store the index of the point in allPoints rather than the actual point.
        std::vector<std::vector<int> > clusters;

	//stores distances of samples from centers
	double** distance; 

	//stores the standard deviation vector for each sample subset (each column of
	//this array corresponds to a particular cluster's standard deviation vector)
	double** STDVector;

	//stores the maximum element of each column in STDVector.
	double* Vmax;

	// stores the index of the component of each column of STDVector which had 
	// the highest value
	int* Vmax_index;

	//stores average distances of samples in each cluster from their
	//corresponding cluster center
	double*  average_distances; 

	//overall average distance.
	double OverallD;

        std::vector<PairDistanceNode> CenterDistances;

        Point* points_helper(int PointCount);
	
	// searches for a point in centers vector;
	int   find_center(Point* to_find);
        

};
#endif
