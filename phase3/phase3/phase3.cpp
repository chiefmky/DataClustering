//Moise Mokoy
//Phase 3

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>

using namespace std;

class Kmeans
{
	//Initialization of vectors
	vector <vector <double>> points, file, centroids;
	vector <double> SE, maxSSE, sumK, SSE, sumDist;

	// Define variables
	int K, N, D, I;
	bool skip;// skip a line of repeated code
	double maxSE, T, maxVal, minVal, indexCH, indexSW, indexD, indexDB;



public:
	Kmeans(int K, int N, int D, int I, double T, vector <vector <double>> file);
	void randomSelection();
	void clusters();
	void newCentroids(int iter);
	void randomPartition();
	void maximin();
	void validationCH(int iter);
	void validationSW();
	void validationD();
	void validationDB();

	double getSSE(int iter);
	double getIndexCH();
	double getIndexSW();
	double getIndexD();
	double getIndexDB();

};

/// <summary>
/// read file 
/// </summary>
/// <param name="F">file</param>
/// <param name="N">number of points</param>
/// <param name="D">number of dimension</param>
/// <returns></returns>
vector <vector <double>> readFile(string F, int& N, int& D)
{
	int n = 0;
	ifstream inFile;
	vector <vector <double>> data;

	inFile.open(F.c_str());// open file
	inFile >> N;//read number of points and assigned to our integer N
	inFile >> D;//read number of dimension and assigned t our integer D
	data.resize(N, vector<double>(D, 0));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			inFile >> data[i][j];
		}
	}
	inFile.close();
	return data;
}


vector <vector <double>> minMaxNormalization(vector <vector <double>> file, int N, int D)
{
	vector <vector <double>> data(N, vector<double>(D));
	vector <double> maxVal(D), minVal(D);
	double newMax = 1.0, newMin = 0.0;



	for (int i = 0; i < D; i++) {
		minVal[i] = file[0][i];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			if (maxVal[j] < file[i][j])
			{
				maxVal[j] = file[i][j];
			}

			if (minVal[j] > file[i][j])
			{
				minVal[j] = file[i][j];
			}

		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			if (abs(maxVal[j] - minVal[j]) > 0)
			{
				data[i][j] = ((file[i][j] - minVal[j]) / (maxVal[j] - minVal[j])) * (newMax - newMin) + newMin;
			}
			else
			{
				data[i][j] = 0.0 * (newMax - newMin) + newMin;
			}

		}
	}

	return data;
}


vector <vector <double>> zScoreNormalization(vector <vector <double>> file, int N, int D)
{
	vector <vector <double>> data(N, vector<double>(D));
	vector <double> meanVal(D), stdVal(D), varVal(D), sumVal1(D), sumVal2(D);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			sumVal1[j] += file[i][j];

		}
	}

	for (int i = 0; i < D; i++) {
		meanVal[i] = sumVal1[i] / N;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			sumVal2[j] += (file[i][j] - meanVal[j]) * (file[i][j] - meanVal[j]);

		}
	}

	for (int i = 0; i < D; i++)
	{
		stdVal[i] = sqrt(sumVal2[i] / N);
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			if (abs(stdVal[j]) > 0)
			{
				data[i][j] = (file[i][j] - meanVal[j]) / stdVal[j];
			}
			else
			{
				data[i][j] = 0.0;
			}

		}
	}

	return data;
}


int main()
{
	string F; //filename
	int I; //max number of iterations
	int R; //number of runs
	int N; // number of points
	int D; //dimentionality of each point
	int iter;
	int Kmin, Kmax; //number of cluster
	double T; //convergent threshold
	ofstream outFile;
	vector <vector <double>> originalData, minMaxData, zScoreData;
	


	//Get the input from the user
	cout << "Enter your the File name: ";
	cin >> F;

	originalData = readFile(F, N, D);
	minMaxData = minMaxNormalization(originalData, N, D);
	R = 100;
	I = 100;
	T = 0.001;
	Kmin = 2;
	Kmax = round(sqrt(N / 2));
	vector <double> maxIndexSW(Kmax), maxIndexCH(Kmax), maxIndexD(Kmax), minIndexDB(Kmax, DBL_MAX);
	





	srand(time(NULL));

	 //Create the output file
	size_t lastindex = F.find_last_of(".");
	string filename = F.substr(0, lastindex);
	outFile.open(F + "_out.txt");

	// append input data to the output file
	outFile << "filename: " << F << ", num of iter:" << I << ", Threshold:" << T << ", num of run: " << R << endl << endl;
	cout << "program running ..." << endl << endl;


	outFile << "Min Max Normalization with Maximin." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		for (int k = Kmin; k <= Kmax; k++)
		{
			Kmeans kmeansSW = Kmeans(k, N, D, I, T, minMaxData);
			iter = 0;	

			kmeansSW.maximin();
			kmeansSW.clusters();
			kmeansSW.newCentroids(iter);
			do {
				iter++;
				kmeansSW.clusters();
				kmeansSW.newCentroids(iter);
				
			} while (!((iter == I) || (((kmeansSW.getSSE(iter - 1) - kmeansSW.getSSE(iter)) / kmeansSW.getSSE(iter - 1)) < T)));

			kmeansSW.validationSW();
			if (kmeansSW.getIndexSW() > maxIndexSW[k - 2])
			{
				maxIndexSW[k - 2] = kmeansSW.getIndexSW();
			}
			
		}
	}

	for (int i = 1; i <= R; i++)
	{
		for (int k = Kmin; k <= Kmax; k++)
		{
			Kmeans kmeansCH = Kmeans(k, N, D, I, T, minMaxData);
			iter = 0;
			

			kmeansCH.maximin();
			kmeansCH.clusters();
			kmeansCH.newCentroids(iter);
			do {
				iter++;
				kmeansCH.clusters();
				kmeansCH.newCentroids(iter);	
			} while (!((iter == I) || (((kmeansCH.getSSE(iter - 1) - kmeansCH.getSSE(iter)) / kmeansCH.getSSE(iter - 1)) < T)));

			kmeansCH.validationCH(iter);
			if (kmeansCH.getIndexCH() > maxIndexCH[k - 2])
			{
				maxIndexCH[k - 2] = kmeansCH.getIndexCH();
			}

		}
	}

	for (int i = 1; i <= R; i++)
	{
		for (int k = Kmin; k <= Kmax; k++)
		{
			Kmeans kmeansDB = Kmeans(k, N, D, I, T, minMaxData);
			iter = 0;


			kmeansDB.maximin();
			kmeansDB.clusters();
			kmeansDB.newCentroids(iter);
			do {
				iter++;
				kmeansDB.clusters();
				kmeansDB.newCentroids(iter);
			} while (!((iter == I) || (((kmeansDB.getSSE(iter - 1) - kmeansDB.getSSE(iter)) / kmeansDB.getSSE(iter - 1)) < T)));

			kmeansDB.validationDB();
			
			if (kmeansDB.getIndexDB() < minIndexDB[k - 2])
			{
				minIndexDB[k - 2] = kmeansDB.getIndexDB();
			}

		}
		
	}

	//for (int i = 1; i <= R; i++)
	//{
	//	for (int k = Kmin; k <= Kmax; k++)
	//	{
	//		Kmeans kmeansD = Kmeans(k, N, D, I, T, minMaxData);
	//		iter = 0;


	//		kmeansD.maximin();
	//		kmeansD.clusters();
	//		kmeansD.newCentroids(iter);
	//		do {
	//			iter++;
	//			kmeansD.clusters();
	//			kmeansD.newCentroids(iter);
	//		} while (!((iter == I) || (((kmeansD.getSSE(iter - 1) - kmeansD.getSSE(iter)) / kmeansD.getSSE(iter - 1)) < T)));

	//		kmeansD.validationD();

	//		if (kmeansD.getIndexD() < maxIndexD[k - 2])
	//		{
	//			maxIndexD[k - 2] = kmeansD.getIndexD();
	//		}

	//	}

	//}

	cout << "Best validation index for eack run for all the K value" << endl << endl;

	outFile << "Silhouette Width index" << endl;
	for (int k = 0; k <= Kmax-2; k++)
	{
		outFile << "K=" << k + 2 << " " << maxIndexSW[k] << endl;
	}
	outFile << endl << endl;

	outFile << "Calinski–Harabasz index" << endl;
	for (int k = 0; k <= Kmax - 2; k++)
	{
		outFile << "K=" << k + 2 << " " << maxIndexCH[k] << endl;
	}
	outFile << endl << endl;
	outFile << "the Davies–Bouldin index" << endl;
	for (int k = 0; k <= Kmax - 2; k++)
	{
		outFile << "K=" << k + 2 << " " << minIndexDB[k] << endl;
	}
	/*outFile << endl << endl;
	outFile << "the Dunn index" << endl;
	for (int k = 0; k <= Kmax - 2; k++)
	{
		outFile << "K=" << k + 2 << " " << maxIndexD[k] << endl;
	}*/

	outFile.close();

	return 0;
}


Kmeans::Kmeans(int K, int N, int D, int I, double T, vector <vector <double>> file)
{
	this->K = K;
	this->N = N;
	this->D = D;
	this->I = I;
	this->T = T;
	this->file = file;

	centroids.resize(K, vector<double>(D, 0));
	SE.resize(K, 0);
	sumDist.resize(K, 0);
	SSE.resize(I, 0);
	maxSSE.resize(D, 0);
	sumK.resize(N, 0);
	points.resize(N, vector<double>(K, 0));
}

/// <summary>
/// function for getting initial randomly the centroids
/// </summary>
void Kmeans::randomSelection()
{
	int randNum;
	for (int i = 0; i < K; i++)
	{
		randNum = rand() % N;
		centroids[i] = file[randNum];
	}
}
/// <summary>
///function for getting initial the centroids using Random Partition
/// </summary>
void Kmeans::randomPartition()
{
	double sum;
	int count, randNum;

	randomSelection();//randomly inital cluster

	for (int i = 0; i < N; i++)
	{
		randNum = rand() % K;
		points[i][randNum] = 1;
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			sum = 0;
			count = 0;
			for (int m = 0; m < N; m++)
			{
				if (points[m][i] == 1)
				{
					sum = sum + file[m][j];
					count++;
				}
			}
			centroids[i][j] = sum / count;
		}
	}

}

void Kmeans::maximin()
{
	double sumOfDist, gSE;
	vector <double> min;
	int randNum;

	randNum = rand() % N;

	centroids[0] = file[randNum];

	for (int i = 0; i < N; i++)
	{
		points[i][0] = 1;

		for (int j = 1; j < K; j++)
		{
			points[i][j] = 0;
		}
	}
	skip = true;



	min.resize(N, 0);

	gSE = 0;
	min = file[0];

	for (int i = 1; i < K; i++)
	{
		sumOfDist = 0;

		for (int j = 0; j < N; j++)
		{
			for (int m = 0; m < D; m++)
			{
				if (points[j][i] == 1)
				{
					sumOfDist += ((file[j][m] - centroids[i][m]) * (file[j][m] - centroids[i][m]));
				}
				if (sumOfDist > gSE)
				{
					min = file[j];
					gSE = sumOfDist;
				}
			}
		}

		centroids[i] = min;

	}


}


/// <summary>
/// cluster points to the nearest centroid
/// </summary>
void Kmeans::clusters()
{
	double sumOfDist;
	double minDist;
	vector <vector <double>> dist(N, vector <double>(K));

	//skip this line of code the first time because it's repeated from maximin
	if (skip == false)
	{
		for (int i = 0; i < N; i++)
		{
			points[i][0] = 1;

			for (int j = 1; j < K; j++)
			{
				points[i][j] = 0;
			}
		}
	}

	skip = false;
	

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < K; j++)
		{
			sumOfDist = 0;
			for (int m = 0; m < D; m++)
			{
				sumOfDist += ((file[i][m] - centroids[j][m]) * (file[i][m] - centroids[j][m]));
			}
			dist[i][j] = sumOfDist;
		}
	}

	for (int i = 0; i < N; i++)
	{
		minDist = dist[i][0];
		for (int j = 1; j < K; j++)
		{
			if (dist[i][j] < minDist)
			{
				minDist = dist[i][j];
				for (int m = 0; m < K; m++)
				{
					points[i][m] = 0;
				}
				points[i][j] = 1;
			}
		}
	}

	for (int i = 0; i < K; i++)
	{
		sumK[i] = 0;
		for (int j = 0; j < N; j++)
		{
			sumK[i] += points[j][i];
		}
	}
}

/// <summary>
/// get new centroids
/// </summary>
/// <param name="iter"></param>
void Kmeans::newCentroids(int iter)
{
	double sum, sumOfDist;
	int count;

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			sum = 0;
			count = 0;
			for (int m = 0; m < N; m++)
			{
				if (points[m][i] == 1)
				{
					sum = sum + file[m][j];
					count++;
				}
			}
			centroids[i][j] = sum / count;
		}
		// check for empty centroid
		if (sumK[i] == 0)
		{
			centroids[i] = maxSSE;
		}

	}

	maxSE = 0;
	maxSSE = file[0];

	for (int i = 0; i < K; i++)
	{
		sumOfDist = 0;
		SE[i] = 0;

		for (int j = 0; j < N; j++)
		{
			for (int m = 0; m < D; m++)
			{
				if (points[j][i] == 1)
				{
					sumOfDist += (file[j][m] - centroids[i][m]) * (file[j][m] - centroids[i][m]);
				}
				if (sumOfDist > maxSE)
				{
					maxSSE = file[j];
					maxSE = sumOfDist;
				}
			}
		}
		SE[i] = sumOfDist;
		SSE[iter] += SE[i];
	}
}

void Kmeans::validationCH(int iter)
{
	vector <vector <double>> meanC(K, vector<double>(D)); //means cluster
	vector <vector <double>> meanT(D, vector<double>(K)); //mean cluster transpose
	vector <double> mean(D);
	vector <int> count(K);
	double trSB = 0;
	double trSW = SSE[iter];


	//claculate mean for all points
	for (int i = 0; i < D; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mean[i] += file[j][i];
		}
		mean[i] = mean[i] / N;
	}

	// mean of cluster point i
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			count[i] = 0;
			for (int m = 0; m < N; m++)
			{
				if (points[m][i] == 1)
				{
					meanC[i][j] += file[m][j];
					count[i]++;
				}
			}
			meanC[i][j] = (meanC[i][j] / count[i]) - mean[j];
			meanT[j][i] = meanC[i][j];//Transpose
		}
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			if (i == j)//only change diagonal value
			{
				trSB += count[i] * meanC[i][j] * meanT[j][i];
			}
		}
	}

	indexCH = ((N - K) / (K - 1)) * (trSB / trSW);
}

void Kmeans::validationSW()
{
	vector <double> meanMinOut(N), meanIn(N);
	vector <vector <double>> dist(N, vector<double>(K, 0));
	double sumOfDist, minDist, maxMean;
	int count = 0;
	double siCoef = 0; //Silhouette Coefficient

	for (int i = 0; i < N; i++)
	{
		if (points[i][0] == 0)
		{
			dist[i][0] = 1;
			sumOfDist = 0;
			for (int j = 0; j < D; j++)
			{
				sumOfDist += (file[i][j] - centroids[0][j]) * (file[i][j] - centroids[0][j]);
			}
			minDist = sumOfDist;
		}
		else
		{
			dist[i][1] = 1;
			sumOfDist = 0;
			for (int j = 0; j < D; j++)
			{
				sumOfDist += (file[i][j] - centroids[1][j]) * (file[i][j] - centroids[1][j]);
			}
			minDist = sumOfDist;
		}
		for (int j = 0; j < K; j++)
		{
			sumOfDist = 0;
			if (points[i][j] != 1)
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += (file[i][m] - centroids[j][m]) * (file[i][m] - centroids[j][m]);
				}
				if (sumOfDist < minDist)
				{
					minDist = sumOfDist;
					for (int n = 0; n < K; n++)
						dist[i][n] = 0;
					dist[i][j] = 1;
				}
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		sumOfDist = 0;
		for (int j = 0; j < N; j++)
		{
			if (i != j && points[j] == dist[i])
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += abs(file[i][m] - file[j][m]);
					count++;
				}
			}
		}
		meanMinOut[i] = sumOfDist / count;
	}
	count = 0;
	for (int i = 0; i < N; i++)
	{
		sumOfDist = 0;
		for (int j = 0; j < N; j++)
		{
			if (i != j && points[i] == points[j])
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += abs(file[i][m] - file[j][m]);
					count++;
				}
			}
		}
		meanIn[i] = sumOfDist / (count - 1);
	}

	for (int i = 0; i < N; i++)
	{
		if (meanMinOut[i] > meanIn[i])
			maxMean = meanMinOut[i];
		else
			maxMean = meanIn[i];
		siCoef += (meanMinOut[i] - meanIn[i]) / maxMean;
	}

	indexSW = siCoef / N;
}

void Kmeans::validationD()
{
	vector <vector <double>> dist(N, vector<double>(K, 0));
	double sumOfDist, minDist, maxDist;
	int count = 0;

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			sumOfDist = 0;
			if (i != j)
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += centroids[i][m] - centroids[j][m];
				}
				minDist = sumOfDist;
			}
			
		}
	}


}
void Kmeans::validationDB()
{
	vector <int> count(K);
	double DB = 0;
	vector <double> maxDBij(K), clusterMean(K), dispersionVal(K);
	vector <vector <double>> dispersion(K, vector<double>(D)), meanC(K, vector<double>(D)), DBij(K, vector<double>(K));


	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			count[i] = 0;
			for (int m = 0; m < N; m++)
			{
				if (points[m][i] == 1)
				{
					meanC[i][j] += file[m][j];
					count[i]++;
				}
			}
			meanC[i][j] = (meanC[i][j] / count[i]);
			clusterMean[i] += (meanC[i][j]);
		}
		clusterMean[i] = clusterMean[i] / D;
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			count[i] = 0;
			for (int m = 0; m < N; m++)
			{
				dispersion[i][j] += abs((file[m][j] - meanC[i][j]) * (file[m][j] - meanC[i][j]));
				count[i]++;
			}
			dispersion[i][j] = sqrt(dispersion[i][j] / count[i]);
			dispersionVal[i] += dispersion[i][j];
		}
		dispersionVal[i] = dispersionVal[i] / D;
	}


	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (j != i)
			{
				if (clusterMean[i] != clusterMean[j])
					DBij[i][j] = (dispersionVal[i] + dispersionVal[j]) / abs(clusterMean[i] - clusterMean[j]);
			}
		}
	}

	maxDBij[0] = DBij[0][1];

	for (int i = 0; i < K; i++)
	{
		for (int j = 1; j < K; j++)
		{
			maxDBij[i] = DBij[i][0];
		}
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (i != j && DBij[i][j] > maxDBij[i])
			{
				maxDBij[i] = DBij[i][j];
			}
		}
		DB += maxDBij[i];
	}

	indexDB = DB / K;
}

//get SSE value
double Kmeans::getSSE(int iter) {
	double val = SSE[iter];
	return val;
}

double Kmeans::getIndexCH() {
	return indexCH;
}

double Kmeans::getIndexSW() {
	return indexSW;
}

double Kmeans::getIndexD() {
	return indexD;
}

double Kmeans::getIndexDB() {
	return indexDB;
}

