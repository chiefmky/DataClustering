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
	double maxSE, T, maxVal, minVal;
	


public:
	Kmeans(int K, int N, int D, int I, double T, vector <vector <double>> file);
	void randomSelection();
	void clusters();
	void newCentroids(int iter);
	void randomPartition();
	void maximin();
	double getSSE(int iter);
	
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
	int K; //number of cluster
	int I; //max number of iterations
	int R; //number of runs
	int N; // number of points
	int D; //dimentionality of each point
	int iter;
	double T; //convergent threshold
	ofstream outFile;
	vector <vector <double>> originalData, minMaxData, zScoreData;

	//Get the input from the user
	cout << "Enter your inputs: ";
	cin >> F >> K >> I >> T >> R;
	
	originalData = readFile(F, N, D);
	minMaxData = minMaxNormalization(originalData, N, D);
	zScoreData = zScoreNormalization(originalData, N, D);


	srand(time(NULL));

	// Create the output file
	size_t lastindex = F.find_last_of(".");
	string filename = F.substr(0, lastindex);
	outFile.open(F + "_out.txt");

	// append input data to the output file
	outFile << F << " " << K << " " << I << " " << T << " " << R << endl << endl;
	cout << "program running ..." << endl << endl;

	/// <summary>
	/// Iniialize Centroid using Random selection
	/// </summary>
	/// <returns></returns>
	outFile << "No Normalization with Random Selection." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, originalData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomSelection();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "Min Max Normalization with Random Selection." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, minMaxData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomSelection();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "z-score Normalization with Random Selection." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, zScoreData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomSelection();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	/// <summary>
	/// Initialize centroids using Random Partition
	/// </summary>
	/// <returns></returns>

	outFile << "No Normalization with Random Partition." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, originalData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomPartition();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "Min Max Normalization with Random Partition." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, minMaxData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomPartition();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "z-score Normalization with Random Partition." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, zScoreData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.randomPartition();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

/// <summary>
/// Initialize centroids using Maximin
/// </summary>
/// <returns></returns>

	outFile << "No Normalization with Maximin." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, originalData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.maximin();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "Min Max Normalization with Maximin." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, minMaxData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.maximin();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

	outFile << "z-score Normalization with Maximin." << endl << endl;
	for (int i = 1; i <= R; i++)
	{
		Kmeans kmeans = Kmeans(K, N, D, I, T, zScoreData);
		iter = 0;

		outFile << "Run " << i << endl;
		outFile << "-----" << endl;
		kmeans.maximin();
		kmeans.clusters();
		kmeans.newCentroids(iter);
		outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		do {
			iter++;
			kmeans.clusters();
			kmeans.newCentroids(iter);
			outFile << "Iteration " << iter + 1 << ": " << kmeans.getSSE(iter) << endl;
		} while (!((iter == I) || (((kmeans.getSSE(iter - 1) - kmeans.getSSE(iter)) / kmeans.getSSE(iter - 1)) < T)));

		outFile << endl;
	}

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

	//dist.resize(N, vector <double>(K, 0));

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

/// <summary>
/// get Sum of square error value
/// </summary>
/// <param name="iter"></param>
/// <returns></returns>
double Kmeans::getSSE(int iter) {
	double val = SSE[iter];
	return val;
}

