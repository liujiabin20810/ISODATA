//=== File Prolog ================================================================
//      This code was developed by NASA, Goddard Space Flight Center, Code 588
//--- Contents -------------------------------------------------------------------
//      File 'Image.cc' contains the implementaion of 'Image' class as described in
//      'Image.h' in C++.
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
//   May 24, 2002     Nargess Memarsadeghi  NASA GSFC, Code 588
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
//
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
//=== End File Prolog=======================================================

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>  //for 'sort' function call

#include "Image.h"
#include "Point.h"
//#include "KMrand.h"

#define K 0.5
#define  MINCLUS 1
#define  MAXCLUS 255

extern int SQUARED;

/**********************************************************************/
/* Constructor                                                        */
/**********************************************************************/
Image::Image()
{

	io_image = NULL;
	distance = NULL;
	STDVector = NULL;
	Vmax = NULL;
	Vmax_index = NULL;
	ImageSizeInByte = 0;
	NumBands = 0;
	NumClusters = 0;
	NumSamples = 0;
	Deleted = false;
}
/**********************************************************************/
/* Constructor							      */
/**********************************************************************/
Image::Image(int row, int col, int bands, int NumClus, int SAMPRM)
{
	int i;

	NumRow = row;
	NumCol = col;
	NumBands = bands;
	NumClusters = NumClus;
	NumSamples = SAMPRM;
	Deleted = false;
	ImageSizeInByte = row*col;
	STDVector = NULL;
	Vmax = NULL;
	Vmax_index = NULL;
	io_image = NULL;


	if (!is_rand)
	{
		io_image = new unsigned char*[bands];
		if (!io_image)
		{
			std::cerr << "Error 1:Memory Allocation Faield in Image::Image" << std::endl;
			exit(1);
		}

		for (i = 0; i < bands; i++)
		{
			io_image[i] = new unsigned char[ImageSizeInByte];
			if (!io_image[i])
			{
				std::cerr << "Error 2: Memory Allocation Faield in Image::Image" << std::endl;
				exit(1);
			}

		}

		/* Initialize */
		for (i = 0; i < bands; i++)
			for (int j = 0; j < ImageSizeInByte; j++)
				io_image[i][j] = 0;
	}

	allPoints = new Point*[ImageSizeInByte];
	distance = new double*[ImageSizeInByte];
	if (!allPoints || !distance)
	{
		std::cerr << "Error 3: Memory Allocation Faield in Image::Image" << std::endl;
		exit(1);
	}


	for (i = 0; i < ImageSizeInByte; i++)
	{
		allPoints[i] = NULL;
		distance[i] = NULL;
	}


}
/******************************************************************/
/* Destructor							  */
/*****************************************************************/
Image::~Image()
{
	int i, num = centers.size();


	if (io_image != NULL)
	{
		for (i = 0; i < NumBands; i++)
			delete[] io_image[i];


		delete[] io_image;
	}



	for (i = 0; i < ImageSizeInByte && !is_rand; i++)
	{
		delete allPoints[i];

	}

	for (i = 0; i < ImageSizeInByte; i++)
	{
		if (distance != NULL)
			delete[] distance[i];

	}


	for (i = 0; i < num; i++)
		delete centers[i];



	centers.clear();
	delete[] allPoints;
	delete[] distance;
	delete[] average_distances;
	delete[] STDVector;

	io_image = NULL;
	allPoints = NULL;
	distance = NULL;

}
/**********************************************************************/
/* If points have been generated already, just set them to all points */
/**********************************************************************/
void Image::setPoints(KMpointArray all)
{

	for (int i = 0; i < ImageSizeInByte; i++)
	{
		allPoints[i] = new Point(NumBands, all[i]);
		//allPoints[i]->setPoint(all[i]);

	}

}
/*****************************************************************/
/* readImage:							 */
/* This file will read each bands information from its file, and */
/* stores it in a row of the 'image' member			 */
/*****************************************************************/
void Image::readImage(std::string* image_names)
{
	int i;

	for (i = 0; i < NumBands; i++)
	{

		std::ifstream OneBand(image_names[i].c_str(), std::ios::in | std::ios::binary);
		if (!OneBand)
		{
			std::cerr << "File ' " << image_names[i] << " ' could not be opened." << std::endl;
			exit(1);
		}

		OneBand.seekg(0, std::ios::beg);
		OneBand.read((char*)io_image[i], ImageSizeInByte);
		//size=OneBand.tellg();
		OneBand.close();
	}

}
/**********************************************************************/
/* outputs the classified data in one row of the 'io_image' array     */
/* and outputs the result into seperate file.                         */
/**********************************************************************/
void Image::writeClassifiedImage(std::string output_name)
{
	int i, j;


	int num_clusters = clusters.size();
	int size;

	for (i = 0; i < num_clusters; i++)
	{
		size = clusters[i].size();
		for (j = 0; j < size; j++)
			io_image[0][clusters[i][j]] = i + 1;

	}

	std::ofstream OneBand(output_name.c_str(), std::ios::out);
	OneBand.write((char*)io_image[0], ImageSizeInByte);
	OneBand.close();

}
/********************************************************************************/
/* This function looks at the io_image (which its dimensions are number of      */
/* points in image times NumBands), and returns a Point object whose coordinate */
/* are that of the 'PointCount' column in io_image.                             */
/********************************************************************************/
Point* Image::points_helper(int PointCount)
{
	Point* to_return;
	double * p = new double[NumBands];
	if (!p)
	{
		std::cout << "Error 4: Memory Allocation Failed" << std::endl;
		exit(1);

	}
	for (int i = 0; i < NumBands; i++)
	{
		p[i] = io_image[i][PointCount];

	}

	to_return = new Point(NumBands, p);
	if (!to_return)
	{
		std::cout << "Error 5: Memory Allocation Failed" << std::endl;
		exit(1);
	}
	return to_return;

}
/*******************************************************************/
Point* Image::getPoint(int row, int col)
{

	int  index = row*NumCol + col;

	if (index > -1 && index <ImageSizeInByte)
		return allPoints[index];

	return NULL;

}
/*******************************************************************/
Point* Image::getPoint(int index)
{

	if (index>-1 && index < ImageSizeInByte)
		return allPoints[index];

	return NULL;

}
/*******************************************************************/
int Image::size()
{
	return ImageSizeInByte;

}
/*******************************************************************************/
/*calculate the distances of all points from centers, and look up              */
/* these distances later 						       */
/******************************************************************************/
void Image::CalculateDistances()
{
	int i, j;
	unsigned int m;

	int size = centers.size();


	for (i = 0; i < ImageSizeInByte; i++)
	{
		if (Deleted && distance[i] != NULL)
			delete[] distance[i];
		distance[i] = new double[size];
		if (!distance[i])
		{
			std::cout << "Error 6: Memory Allocation Failed" << std::endl;
			exit(1);
		}
	}

	if (SQUARED)
	{
		for (m = 0; m < samples.size(); m++)
		{
			i = samples[m];

			for (j = 0; j < size; j++)
				distance[i][j] = allPoints[i]->Norm2DistanceSquared(centers[j]);
		}
	}
	else
	{
		for (m = 0; m < samples.size(); m++)
		{
			i = samples[m];
			for (j = 0; j < size; j++)
				distance[i][j] = sqrt(allPoints[i]->Norm2DistanceSquared(centers[j]));
		}

	}

}
/******************************************************************************/
/* sampleCenters selects a set of initial centers from the set of allPoints.     */
/******************************************************************************/
void Image::sampleCenters()
{
	srand(sample_seed);
	int i, to_add;
	for (i = 0; i < NumClusters; i++)
	{
		to_add = rand() % ImageSizeInByte;

		//If the random point is not a duplicate, add it to the list of centers.
		Point* p;
		if (is_rand)
			p = new Point(*allPoints[to_add]);
		else
			p = points_helper(to_add);
		if (!p)
		{
			std::cout << "Error 7: Memory Allocation Failed" << std::endl;
			exit(1);
		}

		if (find_center(p) == -1)
		{
			centers.push_back(p);
		}
		else
			i--;

	}//for


}
/*******************************************************************************/
/* here we sample points randomly to perform the iterative clustering on them  */
/* for more information on the method used please see:                         */
/* (1) "Programming Pearls", Addison Wesley, 1986,                             */
/* (2) "More Programming Pearls", Addison Wesley, 1988.                        */
/****************************************************************************/

void Image::samplePoints(double s)
{
	srand(sample_seed);
	int n = ImageSizeInByte;

	double probability = s / n;

	for (int i = 0; (i < ImageSizeInByte) && (s > 0); i++)
	{

		samples.push_back(i);
		//generate a random number between 0 and 1
		// 		double x = ((double)rand() / (double)(RAND_MAX + 1));
		// 
		// 		if (x <= probability)
		// 		{
		// 			samples.push_back(i);
		// 			if (!is_rand)
		// 				allPoints[i] = points_helper(i);
		// 			s -= 1;
		// 		}
		// 
		// 		n = n - 1;
		// 		probability = s / n;
	}

	std::cout << "sample size: " << samples.size() << std::endl;
}
/*******************************************************************************/
/* we just modify the data points for real data by adding a small number in    */
/* order of 0.001 such that points that lie on cell borders of the kd-tree be  */
/* assigned to on random to its children's cell.                               */
/*******************************************************************************/
void Image::addNoise()
{
	int s = samples.size();

	for (int i = 0; i < s; i++)
		*(allPoints[samples[i]]) = *(allPoints[samples[i]])/*+ 0.001*kmRanUnif(-1,1)*/;

}
/*******************************************************************************/
bool Image::WasDeleted()
{
	return Deleted;
}
/***********************************************************************************/
/* Given all the points in the image we take sampled points to determine	   */
/* cluster means, for having the final classified image, one can run this function  */
/* with increment=1 to have all the points classified				   */
/***********************************************************************************/
void Image::PutInCluster()
{
	int i, j, size = centers.size();
	int min_index;
	unsigned int m;

	//reset Deleted.
	Deleted = false;

	// initializing the vector of vectors (vector of clusters)

	clusters.erase(clusters.begin(), clusters.end());
	for (i = 0; i < size; i++)
	{
		std::vector<int> v;
		clusters.push_back(v);

	}

	for (m = 0; m < samples.size(); m++)
	{
		i = samples[m];
		min_index = 0;
		for (j = 1; j < size; j++)
		{
			if (distance[i][j] < distance[i][min_index])
			{
				min_index = j;
			}

		}//for loop j

		//add the i-th point to the cluster whose center was closest to the i-th point.   
		clusters[min_index].push_back(i);

	}//for loop i

#if 0
	for (i=0 ; i < clusters.size(); i++)
	{
		std::cout<<"cluster "<<i<<"size : "<<clusters[i].size()<<" center ";
		centers[i]->print(&cout);

	}
#endif
}


/*******************************************************************************/
/* searches in the 'centers' vector to see if it can find a particular point   */
/* It returns the index of the point in the vector if found, and -1 otherwise  */
/*******************************************************************************/
int Image::find_center(Point* to_find)
{
	int value = -1;
	int size = centers.size(), i;

	for (i = 0; i < size; i++)
	{
		if ((*centers[i]) == (*to_find))
		{
			value = i;
			break;
		}

	}

	return value;
}
/************************************************************************************/
/* This function checks number of points in each cluster. If number of points in    */
/* any cluster is less than NumSamples (desired minimum number of points in each    */
/* cluster), then that cluster and its center gets deleted (Note: only the cluster  */
/* from 'clusters' vector gets deleted, the points in the cluster do not get deleted*/
/* from set of sampled points for further iterative clustering.			    */
/************************************************************************************/
void Image::PostAnalyzeClusters()
{
	int i, count = 0;

	std::vector<Point*>::iterator centers_it = centers.begin();
	std::vector< std::vector <int> >::iterator clusters_it = clusters.begin();

	// If a cluster's size is less than desired, delete that cluster and center.
	for (i = 0; i < (int)clusters.size(); i++)
	{
		if (clusters_it->size() < NumSamples)
		{
			clusters_it = clusters.erase(clusters_it);
			delete centers[i];
			centers_it = centers.erase(centers_it);
			i--;  count++;
			Deleted = true;
		}
		else
		{
			clusters_it++;
			centers_it++;
		}
	}
	if (Deleted)
		std::cout << "\tDeleted " << count << " cluster(s)." << std::endl;


}
/****************************************************************************************/
/* Update each cluster center by setting it to the sample mean of its corresponding set */
/****************************************************************************************/
void Image::UpdateCenters()
{
	int i, num_clusters = clusters.size();
	Point c;


	//centers.erase(centers.begin(), centers.end());

	for (i = 0; i < num_clusters; i++)
	{
		int size = clusters[i].size();
		if (size != 0)
		{
			c = (*allPoints[clusters[i][0]]);
		}
		// add all the points in a cluster.	
		for (int j = 1; j < size; j++)
		{
			c = c + (*allPoints[clusters[i][j]]);

		}
		// find average of all points in a cluster.
		if (size != 0)
		{
			c = c / size;
			if (find_center(&c) == -1)
			{
				(*centers[i]) = c;

			}
		}
	}


}
/************************************************************************************/
/* we need to calculate the average distance of all points in a cluster from that   */
/* cluster center.  In order to do so, we need to look up the right distances from  */
/* distances array. Note that if the runtime parameter SQUARED==1, then 	    */
/* 'average_distances' will contain average of squared distances.		    */
/***********************************************************************************/
void Image::CalculateAverageDistances()
{
	int size = clusters.size();
	int i, j;
	double sum;

	average_distances = new double[size];
	if (!average_distances)
	{
		std::cerr << "Memory allocation for 'average_distance' failed." << std::endl;
		std::cerr << "Exiting the program..." << std::endl;
		exit(1);
	}

	for (i = 0; i < size; i++)
	{
		sum = 0;
		int count = clusters[i].size();
		for (j = 0; j < count; j++)
		{
			//need to add the distances of point at position clusters[i][j] in allPoints from
			//center i (since we are adding distances for cluster i).  Note that the structure of
			//distance array is as follows: 
			//distances[point_num][center_num]: 
			//the distance of the point at position 'point_num' in 'allPoints' array from the
			//point at position 'center_num' in 'centers' vector.

			sum += distance[clusters[i][j]][i];
		}
		if (count != 0)
		{
			// NM: Just added 'sqrt' below to make it match with eff2 version
			if (SQUARED)
				average_distances[i] = sqrt(sum / count);
			else
				average_distances[i] = sum / count;
		}
		else
			average_distances[i] = -1;
	}

}
/************************************************************************************/
/* This function retunrs:							    */
/* a) over all average distances of all points in the sample set form their closest */
/*    center if SQUARED==0							    */
/* b) over all average of squared distances of all points in the sample set form    */
/*    their closest center if SQUARED==1					    */
/************************************************************************************/
double Image::OverallAverageDistances()
{
	double sum = 0;
	int i, size = clusters.size();

	for (i = 0; i < size; i++)
		sum += ((clusters[i].size())*average_distances[i]);

	OverallD = sum / samples.size();

	return OverallD;
}
/**************************************************************************************/
/* This function calculated average of sum of *squared* distances of all points from  */
/* their cluster centers.  This function gets called as a basis for comaprison of     */
/* results of this algorithm verses other clustering algorithms that use squared      */
/* distances . This, function will be called at the end of the clustering, and so     */
/* modifications to arrays should not affect performance of the clustering.           */
/**************************************************************************************/
double Image::getDistortions()
{
	int size = clusters.size();
	int i, j;
	double sum = 0;
	if (!SQUARED)
	{
		SQUARED = 1;
		CalculateDistances();
	}

	for (i = 0; i < size; i++)
	{
		int count = clusters[i].size();
		for (j = 0; j < count; j++)
		{
			//need to add the distances of point at position clusters[i][j] in allPoints from
			//center i (since we are adding distances for cluster i).  Note that the structure of
			//distance array is as follows:
			//distances[point_num][center_num]:
			//the distance of the point at position 'point_num' in 'allPoints' array from the
			//point at position 'center_num' in 'centers' vector.

			sum += distance[clusters[i][j]][i];
		}
	}


	OverallD = sum / samples.size();

	return OverallD;

}
/*************************************************************************************/
/* given the position of a point in the image array, we find its coordinate, assuming*/
/* the first points coordinate is (1,1).					     */
/*************************************************************************************/
void Image::printCoordinates(int pos)
{
	std::cout << "( ";
	if ((pos % NumCol) != 0)
		std::cout << (pos / NumCol) + 1;
	else
		std::cout << pos / NumCol;

	std::cout << " ,";

	if ((pos % NumCol) != 0)
		std::cout << pos%NumCol;
	else
		std::cout << NumCol;

	std::cout << " )" << std::endl;

}
/***********************************************************************************/
int Image::getNumCenters()
{
	return centers.size();

}
/*****************************************************************************************/
/* This function calculates a vector for each cluster (a point for each cluster).  	 */
/* Each coordinate D of a particular  cluster's standard deviation vector is the 	 */
/* standard deviation of all coordinate D-s of all points in that cluster from coordinate*/
/* D of the same cluster's center.							 */
/*****************************************************************************************/
void Image::CalculateSTDVector()
{

	int num_clusters = clusters.size();
	int i, b, c, cluster_size;
	double sum, avg;

	if (STDVector == NULL)
	{
		STDVector = new double*[num_clusters];

		if (!STDVector)
		{
			std::cout << "Error 8: Memory Allocation Failed" << std::endl;
			exit(1);
		}
		for (i = 0; i < num_clusters; i++)
		{
			STDVector[i] = new double[NumBands];
			if (!STDVector[i])
			{
				std::cout << "Error 9: Memory Allocation Failed" << std::endl;
				exit(1);
			}


		} //for

	} //if


	for (c = 0; c < num_clusters; c++)
	{
		sum = 0; avg = 0;
		cluster_size = clusters[c].size();


		for (b = 0; b < NumBands; b++)
		{
			for (i = 0, sum = 0; i < cluster_size; i++)
			{
				double * p = allPoints[clusters[c][i]]->getPoint();
				// sum+=pow( (io_image[b][clusters[c][i]]-centers[c]->getCoordinate(b)), 2);
				sum += pow((p[b] - centers[c]->getCoordinate(b)), 2);
			}
			if (cluster_size != 0)
				avg = sum / cluster_size;

			STDVector[c][b] = sqrt(avg);

		}// for b
	}// for c


}
/***************************************************************************/
/* calculates the maximum element in each column of STDVector. Since after */
/* this step we do not need values of STDVector, and each IsoClus iteration*/
/* needs to recalculate STDVector, we will delete the allocated memory for */
/* STDVector at this function as well.					   */
/**************************************************************************/
void Image::CalculateVmax()
{
	int c, b;
	int num_clusters = clusters.size();

	if (Vmax != NULL)
	{
		delete[] Vmax;
		Vmax = NULL;

	}

	if (Vmax_index != NULL)
	{
		delete[] Vmax_index;
		Vmax_index = NULL;

	}

	Vmax = new double[num_clusters];
	if (!Vmax)
	{
		std::cout << "Memory Allocation for 'Vmax' Failed." << std::endl;
		std::cout << "Exitting the program..." << std::endl;
		exit(1);
	}

	Vmax_index = new int[num_clusters];
	if (!Vmax_index)
	{
		std::cout << "Memory Allocation for 'Vmax_index' Failed." << std::endl;
		std::cout << "Exitting the program..." << std::endl;
		exit(1);
	}



	for (c = 0; c < num_clusters; c++)
	{
		Vmax[c] = STDVector[c][0];
		Vmax_index[c] = 0;
		for (b = 1; b < NumBands; b++)
		{
			if (STDVector[c][b] > Vmax[c])
			{
				Vmax[c] = STDVector[c][b];
				Vmax_index[c] = b;
			}
		}

		delete[] STDVector[c];

	}
	delete[] STDVector;
	STDVector = NULL;
}
/**************************************************************************************/
/* This function evaluates 3 conditions X, A, and B and evaluates the value of:	      */
/* (X and (A or B)) for each cluster.						      */
/* These conditions are described in step 10 of IsoClus algorithm.  See top of        */
/* IsoClus.cc for references.							      */
/* It adds the cluster numbers of all those clusters that satisfied (X and (A or B))  */
/**************************************************************************************/
std::vector<int> Image::ShouldSplit(double stdv)
{
	std::vector<int> to_return;
	bool X = false, A = false, B = false;
	int num_clusters = clusters.size();

	//evaluateing X
	for (int i = 0; i < num_clusters; i++)
	{
		if (Vmax[i] > stdv)
			X = true;
		//the goal is to evaluate [ X and (A or B) ]
		if (X)
		{
			//evaluate A
			if ((average_distances[i] > OverallD) && (clusters[i].size() > 2 * (NumSamples + 1)))
				A = true;

			if (A)
				to_return.push_back(i);

			else
			{
				//evaluate B: compares the actual number of clusters with the desired number of
				//clusters.	
				if (num_clusters <= (NumClusters / 2))
					B = true;

				if (B)
					to_return.push_back(i);

			}//else

		} //if (X)

		// reset the booleans for next iteration
		X = false; A = false; B = false;
	} //for loop
	return to_return;
}
/****************************************************************************************/
/* It splits all those clusters whose index appears in 'to_split' parameter.  It does   */
/* so by calculating two new centers.  It assigns one of the centers to the same cluster*/
/* that is being split, and add the other new center to the end of 'centers' array.     */
/****************************************************************************************/
void Image::Split(std::vector<int> to_split)
{
	int i, size = to_split.size(), j, index, num_centers = centers.size();
	double Gj, current;
	Point Zminus, Zplus;
	try
	{
		for (j = 0; j < size && num_centers < MAXCLUS; j++, num_centers++)
		{
			index = to_split[j];
			Gj = K*Vmax[index];


			//Zminus=(*centers[index])-Gj;
			//Zplus=(*centers[index])+Gj;

			Zminus = *centers[index];
			Zplus = *centers[index];

			//i would be the coordinate of STDV vector which had the maximum component.
			i = Vmax_index[index];
			current = Zminus.getCoordinate(i);

			Zminus.setCoordinate(i, current - Gj);
			Zplus.setCoordinate(i, current + Gj);



			// this erases the previous center value and updates it by what is called
			// Zplus in in ISOCLUS algorithm step 10, see references in IsoClus.cc

			(*centers[index]) = Zminus;
			centers.push_back(new Point(Zplus));
			std::cout << "\tSplit cluster " << index + 1 << "." << std::endl;

		} // for

	} //try
	catch (std::bad_alloc exception)
	{
		std::cout << "Exception occured in Image::Split function: " << exception.what() << std::endl;
		exit(1);
	}

}
/********************************************************************************************/
void Image::ComputeCenterDistances()
{
	unsigned int i, j, size, count;
	double dist;
	bool Emptied = false;

	size = centers.size();

	if (size*(size-1)/2 != CenterDistances.size())
	{
		CenterDistances.clear();
		Emptied = true;
	}

	if (size == 0 || Emptied)
	{
		for (i = 0; i < size - 1; i++)
			for (j = i + 1; j < size; j++)
			{
				dist = sqrt(centers[i]->Norm2DistanceSquared(centers[j]));
				PairDistanceNode n = { dist, i, j };
				CenterDistances.push_back(n);

			}
	}

	else  //just need to update "dist" values rather than allocating all structure nodes
	{
		//note that 'count' is suppose to be the index number of the next element in
		//CenterDistances array. 
		for (i = 0, count = 0; i < size - 1; i++)
			for (j = i + 1; j < size; j++, count++)
			{
				dist = sqrt(centers[i]->Norm2DistanceSquared(centers[j]));
				CenterDistances[count].dist = dist;

			}
	}

#if 0
	cout<<" in pair wise distances size was "<<size<<endl;
	cout<<"number of iterations was "<<count<<endl;

	for (i=0; i < count; i++)
	{
		std::cout<<i<<": dist= "<<CenterDistances[i].dist<<"  c1: "<<CenterDistances[i].c1<<" c2: "<<CenterDistances[i].c2<< std::endl; 
	}
#endif

}

/************************************************************************************/
/* Searches the list of pair distance nodes and selects those pairs whose distances */
/* from each other are less than parameter LUMP.  Orders the list in ascending order*/
/* and selects the first MAXPAIR pairs as candid pairs for lumping.                 */
/************************************************************************************/
std::vector<PairDistanceNode> Image::FindLumpCandidates(double lump, int MAXPAIR)
{
	int count = 0;
	int size = CenterDistances.size();
	std::vector<PairDistanceNode> lump_candidates;

	// sort the list of center pairs based on '<' operator which is overloaded to compare
	// 'dist' field of each node. (sort based on distances)
	sort(CenterDistances.begin(), CenterDistances.end());

	for (int i = 0; i < size && lump != 0 && count < MAXPAIR; i++)
	{
		if (CenterDistances[i] <lump)
		{
			lump_candidates.push_back(CenterDistances[i]);
			count += 1;
		}
	}

	sort(lump_candidates.begin(), lump_candidates.end());

	//If there are more candidates to lump, than the MAXPAIR, select the first MAXPAIR candidates

	if (count > MAXPAIR)
		lump_candidates.erase(lump_candidates.begin() + MAXPAIR, lump_candidates.end());


	return lump_candidates;
}
/**********************************************************************************/
/* Eligible clusters among to_lump vector will be lumped.                         */
/* Each cluster/center can be lumped only once, and thus not all the centers      */
/* associated with elements of to_lump vector will be lumped.                     */
/**********************************************************************************/
void Image::Lump(const  std::vector<PairDistanceNode>& to_lump)
{
	int to_lump_size = to_lump.size();
	int orig_centers_size = centers.size();
	int i, clus1_size, clus2_size, clus1, clus2, used_index, count = 0;

	try
	{
		bool* used_centers = new bool[orig_centers_size];

		for (i = 0; i < orig_centers_size; i++)
			used_centers[i] = 0;


		for (i = 0; i < to_lump_size; i++)
		{
			clus1 = to_lump[i].c1;
			clus2 = to_lump[i].c2;

			// only go about lumping the two clusters if neither of the centers had been used before.
			if (!used_centers[clus1] && !used_centers[clus2] && (orig_centers_size - count) > MINCLUS)
			{
				//calculate the new center
				clus1_size = clusters[clus1].size();
				clus2_size = clusters[clus2].size();

				Point new_center = ((*centers[clus1]) * clus1_size) + ((*centers[clus2]) * clus2_size);
				new_center = new_center / (clus1_size + clus2_size);
				centers.push_back(new Point(new_center));

				//merge the two clusters into a new one, and add the new cluster to the end of vector of clusters.
				std::vector<int> new_cluster;
				new_cluster.insert(new_cluster.end(), clusters[clus1].begin(), clusters[clus1].end());
				new_cluster.insert(new_cluster.end(), clusters[clus2].begin(), clusters[clus2].end());
				clusters.push_back(new_cluster);

				used_centers[clus1] = 1;
				used_centers[clus2] = 1;
				count++;
				std::cout << "\tLumped clusters " << clus1 + 1 << " and " << clus2 + 1 << "." << std::endl;
			} // if could be lumped


		} // for loop



		std::vector<Point*>::iterator centers_it = centers.begin();
		std::vector<  std::vector <int> >::iterator clusters_it = clusters.begin();

		for (i = 0, used_index = 0; used_index < orig_centers_size; i++, used_index++)
		{
			//if this cluster has been lumped remove it and its center.
			if (used_centers[used_index])
			{
				clusters_it = clusters.erase(clusters_it);
				delete centers[i];
				centers_it = centers.erase(centers_it);
				i--; 
			} 
			else
			{
				centers_it++, clusters_it++;
			}// if
		} // for

		delete[] used_centers;
	} //try
	catch (std::bad_alloc exception)
	{
		std::cout << "Exception occured in Image::Lump function: " << exception.what() << std::endl;
		exit(1);
	}


}
/***************************************************************************************/
/* Select all points for the last iterative clustering.				       */
/***************************************************************************************/
void Image::preFinalClustering()
{
	//since it is going to be the last round, we need to classify every single point
	samples.clear();
	for (int i = 0; i < ImageSizeInByte; i++)
	{
		samples.push_back(i);
		if (allPoints[i] == NULL)
		{
			allPoints[i] = points_helper(i);

		}
	}
}
/***************************************************************************************/
/* Print a report of current clustering results. (Number of clusters, their, sizes, etc*/
/***************************************************************************************/
void Image::generateReport()
{
	int i, num_clusters = clusters.size();

	std::cout << "\n\tNumber of Clusters: " << clusters.size() << std::endl;
	std::cout << "\t==============================================================================\n";
	for (i = 0; i < num_clusters; i++)
	{
		std::cout << "\tCluster " << i + 1 << ":\n\tSize: " << clusters[i].size();
		if (SQUARED)
			std::cout << "\n\tAverage of squared distances: ";
		else
			std::cout << "\n\tAverage of distances: ";
		std::cout << average_distances[i] << std::endl;
		std::cout << "\tCenter:";
		centers[i]->print(&std::cout);
		std::cout << "\t==============================================================================\n";

	}
	if (SQUARED)
		std::cout << "\tOverall average of squared distances of all points from their cluster center: ";
	else
		std::cout << "\tOverall average of distances of all points from their cluster center: ";

	std::cout << OverallD << std::endl;

}

void Image::getClusterLabel(int* label)
{
	if (NULL == label)
		return;

	int i,j, num_clusters = clusters.size();
	int size;
	std::vector<double> dists(num_clusters);
	for (i = 0; i < num_clusters; i++)
	{
		size = clusters[i].size();
		std::vector<double> _dist(size);
		for (j = 0; j < size; j++)
			_dist[j] = distance[clusters[i][j]][i];

		std::sort(_dist.begin(), _dist.end());

		int N = (size-1);
		dists[i] = _dist[N];
//		std::cout << "max distance: "<<i<<" " << dists[i] << std::endl;
	}

	for (i = 0; i < num_clusters; i++)
	{
		size = clusters[i].size();
		for (j = 0; j < size; j++)
		{
			if (distance[clusters[i][j]][i] <= dists[i]) // max distance 0.01
				label[clusters[i][j]] = i + 1;
			else
				label[clusters[i][j]] = 0;
		}
	}
}

