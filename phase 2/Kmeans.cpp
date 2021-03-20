//Moise Mokoy
//Phase 2

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
	vector <double> SE, maxSSE, sumK, SSE;

	// Define variables
	int K, N, D, I;
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
				data[i][j] = ((file[i][j] - minVal[j]) / (maxVal[j] - minVal[j])) * (newMax - newMin) + newMin ;
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
	int I = 100; //max number of iterations
	int R = 100; //number of runs
	int N; // number of points
	int D; //dimentionality of each point
	int iter;
	int Kmin, Kmax; //number of cluster
	double T = 0.001; //convergent threshold
	ofstream outFile;
	vector <vector <double>> originalData, minMaxData, zScoreData;

	//Get the input from the user
	cout << "Enter your the File name: ";
	cin >> F;
	
	originalData = readFile(F, N, D);
	minMaxData = minMaxNormalization(originalData, N, D);
	Kmin = 2;
	Kmax = round(sqrt(N / 2));



	srand(time(NULL));

	// Create the output file
	/*size_t lastindex = F.find_last_of(".");
	string filename = F.substr(0, lastindex);
	outFile.open(F + "_out.txt");*/

	// append input data to the output file
	cout <<"filename: "<< F << ", num of iter:" << I << ", Threshold:" << T << ", num of run: " << R << endl << endl;
	cout << "program running ..." << endl << endl;


	cout << "Min Max Normalization with Maximin." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		for (int k = Kmin; k <= Kmax; k++)
		{
			Kmeans kmeans = Kmeans(K, N, D, I, T, minMaxData);
			iter = 0;

			kmeans.maximin();
			kmeans.clusters();
			kmeans.newCentroids(iter);
			do {
				iter++;
				kmeans.clusters();
				kmeans.newCentroids(iter);
			} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));
			kmeans.validationSW();

			cout << k << " " << kmeans.getIndexSW << endl;

		//	outFile << endl;
		}
			
	}

	//outFile.close();

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
	SE.resize(K,0);
	SSE.resize(I,0);
	maxSSE.resize(D,0);
	sumK.resize(N,0);
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
	double sum, sumOfDist;
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
	double sumOfDist,gSE;
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


	for (int i = 0; i < N; i++)
	{
		points[i][0] = 1;

		for (int j = 1; j < K; j++)
		{
			points[i][j] = 0;
		}
	}

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
					sumOfDist += ((file[j][m] - centroids[i][m]) * (file[j][m] - centroids[i][m]));
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

void Cluster::validationCH(int iter)
{
	vector <vector <double>> meanC(K, vector<double>(D)); //means cluster
	vector <vector <double>> meanT(D, vector<double>(K)); //mean cluster transpose
	vector <double> mean(D);
	vector <int> count(K);
	double trSB = 0;


	for (int i = 0; i < D; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mean[i] += file[j][i];
		}
		mean[i] = mean[i] / N;
	}

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

	//for (int i = 0; i < K; i++)
	//{
	//	for (int j = 0; j < D; j++)
	//	{
	//		meanT[j][i] = meanC[i][j];//Transpose
	//	}
	//}


	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			if (i == j)//only change iagonal value
			{
				trSB += count[z] * meanC[z][j] * meanT[j][z];
			}
		}
	}

	indexCH = ((N - K) * trSB) / ((K - 1) * SSE[iter]);
	//CH = ((N - K) / (K -1)) * (trSB / SSE[iter])
}

void Cluster::validationSW()
{	
	vector <double> meanMinOut(N), meanIn(N);
	vector <vector <double>> dist(N, vector<double>(K, 0));
	double sumOfDist, minDist, max;
	int countIn = 0, countOut = 0;
	double SW = 0;

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
			if (j != i && points[j] == dist[i])
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += abs(file[i][m] - file[j][m]);
					countOut++;
				}
			}
		}
		meanMinOut[i] = sumOfDist / countOut;
	}

	for (int i = 0; i < N; i++)
	{
		sumOfDist = 0;
		for (int j = 0; j < N; j++)
		{
			if (j != i && points[j] == points[i])
			{
				for (int m = 0; m < D; m++)
				{
					sumOfDist += abs(file[i][m] - file[j][m]);
					countIn++;
				}
			}
		}
		meanIn[i] = sumOfDist / (countIn - 1);
	}

	for (int i = 0; i < N; i++)
	{
		if (meanMinOut[i] > meanIn[i])
			max = meanMinOut[i];
		else
			max = meanIn[i];
		SW += ((meanMinOut[i] - meanIn[i]) / max);
	}

	indexSW = SW / N;

}

void Cluster::validationD()
{

}
void Cluster::validationDB()
{
	vector <int> count(K);
	double DB = 0;
	vector <double> maxDBij(K), mu(K), sig(K);
	vector <vector <double>> sigC(K, vector<double>(D)), meanC(K, vector<double>(D)), DBij(K, vector<double>(K));


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
		}
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			mu[i] += (meanC[i][j] * meanC[i][j]);
		}
		mu[i] = mu[i] / D;
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			for (int m = 0; m < N; m++)
			{
				sigC[i][j] += abs((file[m][j] - meanC[i][j]) * (file[m][j] - meanC[i][j]));
			}
			sigC[i][j] = sqrt(sigC[i][j] / count[i]);
		}
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			sig[i] += sigC[i][j];
		}
		sig[i] = sig[i] / D;
	}

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (j != i)
			{
				if (mu[i] != mu[j])
					DBij[i][j] = (sig[i] + sig[j]) / abs(mu[i] - mu[j]);
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
	//cout << DB << endl;
}

/// <summary>
/// get Sum of square error value
/// </summary>
/// <param name="iter"></param>
/// <returns></returns>
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

