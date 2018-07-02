// ISOSample.cpp
// 头文件顺序
//1. 本类声明
//2. C系统文件
//3. C++系统文件
//4. 其他库头文件
//5. 本项目内头文件
#include <vector>
#include <iostream>
#include <random>
#include <map>

#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
// cluster
#include "Image.h"

//log
#include <glog/logging.h>

#ifdef _DEBUG
    #pragma comment(lib,"opencv_core331d.lib")
    #pragma comment(lib,"opencv_highgui331d.lib")
	#pragma comment(lib,"opencv_imgcodecs331d.lib")
    #pragma comment(lib,"opencv_imgproc331d.lib")
	#pragma comment(lib,"glogd.lib")
#else
	#pragma comment(lib,"opencv_core331.lib")
	#pragma comment(lib,"opencv_highgui331.lib")
	#pragma comment(lib,"opencv_imgcodecs331d.lib")
	#pragma comment(lib,"opencv_imgproc331.lib")
	#pragma comment(lib,"glog.lib")
#endif // _DEBUG

int sample_seed = 51;         // for random sampling of centers and points
bool is_rand = true;
int SQUARED = 0;

using  cv::Mat;
using cv::Rect;

KMpointArray AllocPts(int n, int dim)	// allocate n pts in dim
{
	KMpointArray pa = new KMpoint[n];		// allocate points
	/* Nargess: added  memory allocation check */
	if (!pa)
	{
		LOG(ERROR) << "Memory Allocation failed in 'kmAllocPts' for 'pa' ";
		return NULL;
	}

	KMpoint	  p = new KMcoord[n*dim];	// allocate space for coords

	/* Nargess: added  memory allocation check */
	if (!p)
	{
		LOG(ERROR) << "Memory Allocation failed in 'kmAllocPts' for 'p'";
		return NULL;
	}

	for (int i = 0; i < n; i++) {
		pa[i] = &(p[i*dim]);
	}

	return pa;
}

void DeallocPts(KMpointArray &points)
{
	delete[] points[0];				// dealloc coordinate storage
	delete[] points;				// dealloc points
	points = NULL;
}

inline double elapsedTime(clock_t start) {
	return double(clock() - start) / double(CLOCKS_PER_SEC);
}

bool autoCluster(cv::Mat mask)
{
	// prepare data
	// init
//	rects.clear();

	Mat stamp_text_mask;
	mask.copyTo(stamp_text_mask);

	int sampleCount = countNonZero(stamp_text_mask);

	if (sampleCount < 2500)
		return false;

	int w = mask.cols, h = mask.rows;

	int NUMBANDS = 2;
	int SAMPRM = 500;
	int NUMCLUS = 3;

	KMpointArray points = AllocPts(sampleCount, NUMBANDS); // (x,y)

	int ncount = 0;
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			int mask_val = stamp_text_mask.ptr<uchar>(i)[j];
			if (0 == mask_val)
				continue;

			points[ncount][0] = 1.0*j / w;  // x ->j
			points[ncount][1] = 1.0*i / h;  // x->i
			ncount++;
		}

	//////////////////////////////////////////////////////////////////////////
	int iter = 0;
	double exec_time = 0;
	int MAXITER = 20;
	int MAXPAIR = 5;                       // maximum number of pairs to lump
	double LUMP = 0.1;
	double std_m = 1.0 / pow(NUMCLUS, 1 / (double)NUMBANDS);
	double STDV = std_m*0.1;

	Image IMG = Image(sampleCount, 1, NUMBANDS, NUMCLUS, SAMPRM);

	IMG.setPoints(points);
	IMG.sampleCenters();
	IMG.samplePoints(sampleCount);


	clock_t start = clock();                    // start the clock 

	for (iter = 1; iter <= MAXITER; iter++)
	{
		LOG(INFO) << " Iteration Number " << iter << " :";
		if (iter == MAXITER)
		{
			LOG(INFO) << "\tPerform the last iterative clustering on all points";
			IMG.preFinalClustering();
		}

		do
		{
			IMG.CalculateDistances();

			LOG(INFO) << "\tPut points into clusters.";
			IMG.PutInCluster();
			//STEP 3:
			IMG.PostAnalyzeClusters();
			//STEP 4:
			LOG(INFO) << "\tUpdate centers by calculating the average point in each cluster.";
			IMG.UpdateCenters();


		} while (IMG.WasDeleted());

		//need to update distances since in the last iteration centers have modified.

		IMG.CalculateDistances();
		IMG.PutInCluster();
		//STEP 5:
		IMG.CalculateAverageDistances();
		//STEP 6:
		IMG.OverallAverageDistances();
		//STEP 7:
		int next_step = 8;
		if (iter == MAXITER)
		{
			LUMP = 0;
			next_step = 11;
		}
		else if (IMG.getNumCenters() <= (NUMCLUS / 2))
		{
			next_step = 8;
		}
		else if ((iter % 2 == 0) || (IMG.getNumCenters() >= 2 * NUMCLUS))
		{
			next_step = 11;
		}

		switch (next_step)
		{
		case 8:
		{
			//STEP 8:
			IMG.CalculateSTDVector();
			//STEP 9:
			IMG.CalculateVmax();
			//STEP 10:
			// the vector to_split will contain integers that represent the cluster numbers
			// that need to be split.
			std::vector<int> to_split = IMG.ShouldSplit(STDV);
			if (to_split.size() != 0)
			{
				IMG.Split(to_split);
				//we need to substract one if it was the last iteration because otherwise we
				//we will exit the loop without updating clusters.
				if (iter == MAXITER)
					iter = iter - 1;
				break; //go to step 2
			}
		} //CASE 8

		case 11:
		{
			//STEP 11:
			IMG.ComputeCenterDistances();
			//STEP 12:
			std::vector<PairDistanceNode> to_lump = IMG.FindLumpCandidates(LUMP, MAXPAIR);
			//STEP 13:
			if (to_lump.size() != 0)
				IMG.Lump(to_lump);

		} //CASE 11   
		} // SWITCH

		// 		LOG(INFO) << "total overall dist " << IMG.getDistortions();
		if (IMG.getDistortions() < 0.005)
			break;
	}  // for LOOP

	exec_time = elapsedTime(start);  // get elapsed time
	LOG(INFO) << "Algorithm's run time: " << exec_time << " CPU seconds.";
	LOG(INFO) << "total overall dist " << IMG.getDistortions();
	LOG(INFO) << "cluster number " << IMG.getNumCenters();
	//	IMG.generateReport();

	DeallocPts(points);

	int* label = new int[sampleCount];
	if (!label)
	{
		LOG(ERROR) << "Memory Allocation for 'label' Failed.";
		return false;
	}

	int num_cluster = IMG.getNumCenters();
// 	if (num_cluster > 6)
// 	{
// 		LOG(ERROR) << "points cluster more than 6 : num_cluster = " << num_cluster;
// 		return false;
// 	}

	IMG.getClusterLabel(label);
	std::vector< std::vector<cv::Point> > clustered_points(num_cluster);
	ncount = 0;
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			int mask_val = stamp_text_mask.ptr<uchar>(i)[j];
			if (0 == mask_val)
				continue;

			int cluster_id = label[ncount];
			if (0 == cluster_id)
			{
				ncount++;
				continue;
			}
			clustered_points[label[ncount] - 1].push_back(cv::Point(j, i));
			ncount++;
		}

	delete[] label;
	label = NULL;

	//////////////////////////////////////////////////////////////////////////

	cv::Vec3b colors[8] = { cv::Vec3b(0, 0, 255),
		cv::Vec3b(0, 255, 0),
		cv::Vec3b(255, 0, 0),
		cv::Vec3b(255, 0, 255),
		cv::Vec3b(255, 255, 0),
		cv::Vec3b(0, 255, 255),
		cv::Vec3b(255, 255, 255) };

	Mat drawImg;
	cvtColor(mask, drawImg, CV_GRAY2BGR);
	drawImg.setTo(cv::Scalar::all(0));

	std::random_device rd;
	for (int i = 0; i < num_cluster; i++)
	{
		std::mt19937 mt(rd());
		std::uniform_int_distribution<> d(0, 255);

		std::vector<cv::Point> points = clustered_points[i];
		cv::Vec3b color = colors[i%8];

		for (int j = 0; j < points.size(); j++)
		{
			int px = points[j].x, py = points[j].y;
			drawImg.ptr<cv::Vec3b>(py)[px] = color;
		}
	}


	cv::namedWindow("Cluster", 0);
	cv::imshow("Cluster", drawImg);
	cv::waitKey();

	cv::imwrite("cluster.jpg", drawImg);

	return true;
}

void genDataPoints(int width, int height, cv::Mat& image, cv::Point center, int radius)
{

	if (image.empty())
		image = cv::Mat(height,width,CV_8U,cv::Scalar::all(0));

	int total_count = (int)( 0.5+ CV_PI*radius*radius);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<> d(-radius,radius );

	std::map<int, int> hist;
	for (int n = 0; n < total_count*0.8; ++n) {

		int rand1 = d(mt);
		int rand2 = d(mt);

		float distance = sqrt(rand1*rand1 + rand2*rand2);

		int x = center.x + rand1;
		int y = center.y + rand2;

		if (x < 0 || x >= width || y < 0 || y >= height || distance >= radius)
			continue;

		image.ptr<uchar>(y)[x] = 255;

		++hist[rand1];
		++hist[rand2];

	}
	for (auto p : hist) {
		std::cout << std::fixed << std::setprecision(1) << std::setw(2)
			<< p.first << ' ' << std::string(p.second / 20, '*') << '\n';
	}

}
// generate random data.
int main(int argc, char* argv[])
{
	// Init
	google::InitGoogleLogging(argv[0]); 

	FLAGS_logtostderr = true;  // log 打印输出
	FLAGS_log_dir = ""; // log 保存目录
	FLAGS_logbufsecs = 0; // log 实时输出
	FLAGS_colorlogtostderr = true; // 打开log 输出颜色

	LOG(INFO) << "This is a data cluster demo for ISODATA ";
	LOG(INFO) << "=========================";

	int width = 200, height = 200;
	int radius = 50;
	Mat genMask;

	genDataPoints(width, height, genMask,cv::Point(width/2,height/2) ,25);
//	genDataPoints(width, height, genMask, cv::Point(width / 4, height / 4), 25);
//	genDataPoints(width, height, genMask, cv::Point(3*width / 4, 3*height / 4), 25);
	genDataPoints(width, height, genMask, cv::Point(width / 4, 3 * height / 4), 25);
	genDataPoints(width, height, genMask, cv::Point(3*width / 4, height / 4), 25);

	cv::namedWindow("data points", 0);
	cv::imshow("data points", genMask);
	cv::waitKey();
	cv::imwrite("data.jpg", genMask);
	std::vector<cv::Rect> rects;
	autoCluster(genMask);

	// Shut
	google::ShutdownGoogleLogging();
	return 0;
}

