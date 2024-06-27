#include <math.h>
#include <cassert>
#include <string.h>
#include <cerrno>
#include <iostream>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

using namespace std;



/* Print arrays/vectors*/
void printArray(int* input, int columns, int rows);
void printVector(vector<int> input, int columns, int rows);

//Waksman Permutation
int* sampleRandomPermutation(int length);
vector<int> sortingNetworkBits(int* perm, int length);
vector<int> sortingNetworkBits_old(int* perm, int length);
vector<int> WaksmanPermutationSetup(int* perm, int length);
int* computeInversePermutation(int* perm, int length);
int SetSwapper(vector<int> array, int index, int via);
int SetSwapperI(vector<int>* array, int index, int via);
int SetSwapperO(vector<int>* array, int index, int via);
int neighbor(int x);
void evaluateWaksmanNetwork(vector<int> swapbits, int* input, int length);
int* evaluateWaksmanNetwork_old(vector<int> swapbits, int* input, int length);
void swapGate(int* input, int pos1, int pos2, int swapbit);
int count_swapbits(int n); 
