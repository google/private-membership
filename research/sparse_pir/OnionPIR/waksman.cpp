#include "waksman.h"

int count_swapbits(int n){
    int result = 0; 
    for (int i = 1; i <=n; i++){
        result += ceil(log2(i)); 
    }
    return result; 
}

void printArray(int* input, int columns, int rows)
{
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < columns; i++)
		{
			printf("%d", input[i + columns * j]);
			printf(" ");
		}
		cout << "\n";
	}
}

void printVector(vector<int> input, int columns, int rows)
{
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < columns; i++)
		{
			printf("%d", input[i + columns * j]);
			printf(" ");
		}
		cout << "\n";
	}
}




int* sampleRandomPermutation(int length)
{
	srand(time(NULL));
	int passes = length * log2(length);
	int* permutation = (int*)malloc(length * sizeof(int));
	int temp;

	// Initialize permutation to 1 2 3 4 ...
	for (int i = 0; i < length; i++)
	{
		permutation[i] = i;
	}

	//Make enough swaps so that the final vector is close to a random permutation
	for (int k = 0; k < passes; k++)
	{
		for (int i = 0; i < length - 1;i++)
		{
			if ((int)(rand() % 2) == 1)
			{
				temp = permutation[i + 1];
				permutation[i + 1] = permutation[i];
				permutation[i] = temp;
			}
		}
	}

	return permutation;
}

int* computeInversePermutation(int* perm, int length)
{
	int* inverse = (int*)malloc(length * sizeof(int));

	for (int i = 0; i < length; i++)
	{
		inverse[perm[i]] = i;

	}
	return inverse;
}

vector<int> sortingNetworkBits_old(int* perm, int length)
{
	if (length == 1)
	{
		vector<int> swapbits(0);
		return swapbits;
	}
	else if (length == 2)
	{
		vector<int> swapbits(1);
		if (perm[0] == 0)
			swapbits[0] = 0;
		else
			swapbits[0] = 1;
		return swapbits;
	}
	int* inverse = computeInversePermutation(perm, length);
	int ins = (int)length / 2;
	int outs = (int)length / 2 - 1;
	int totalSwaps = length - 1;
	vector<int> I(ins);
	vector<int> O(outs + 1);
	for (int i = 0; i < outs; i++)
		O[i] = -1;
	O[outs ] = 0;
	for (int i = 0; i < ins; i++)
		I[i] = -1;

	int** pi = (int**)malloc(2 * sizeof(int*));
	pi[0] = (int*)malloc(((int)length / 2) * sizeof(int));
	pi[1] = (int*)malloc(((int)length / 2 + (length % 2)) * sizeof(int));

	int j;
	int i;
	int outvia;
	int invia;
	if ((length % 2) == 0)
	{
		j = length - 1;
		i = inverse[j];
		outvia = 1;
		invia = SetSwapperI(&I, i, outvia);
		totalSwaps -= 1;
		pi[outvia][(int)i / 2] = (int)j / 2;
		while (true)
		{
			i = neighbor(i);
			j = perm[i];
			if (j == length - 2)
			{
				pi[0][(int)i / 2] = (int)j / 2;
				break;
			}
			outvia = SetSwapperO(&O, j, invia);
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
			j = neighbor(j);
			i = inverse[j];
			invia = SetSwapperI(&I, i, outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
		}
	}
	else
	{
		j = length - 1;
		i = inverse[j];
		outvia = 1;
		while (i != length - 1)
		{
			invia = SetSwapperI(&I, i, outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
			i = neighbor(i);
			j = perm[i];
			outvia = SetSwapperO(&O, j, invia);
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
			j = neighbor(j);
			i = inverse[j];
		}
		pi[outvia][(int)i / 2] = (int)j / 2;
	}

	int j0 = outs - 1;
	while (totalSwaps > 0) {
		while (O[j0] != -1) {
			j0--;
		}
		j = 2 * j0;
		i = inverse[j];
		while (true) {
			invia = SetSwapperI(&I, i, outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
			i = neighbor(i);
			j = perm[i];
			outvia = SetSwapperO(&O, j, invia);
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
			j = neighbor(j);
			i = inverse[j];
			if (j == 2 * j0)
				break;
		}
	}

	vector<int> pi_up = sortingNetworkBits_old(pi[0], (int)length / 2);
	vector<int> pi_down = sortingNetworkBits_old(pi[1], (int)length / 2 + (length % 2));
	I.insert(I.end(), pi_up.begin(), pi_up.end());
	I.insert(I.end(), pi_down.begin(), pi_down.end());
	I.insert(I.end(), O.begin(), O.end());
	return I;
}

int SetSwapper(vector<int>* swapbits, int index, int via)
{
	if ((index % 2) == 1)
		(*swapbits)[(int)index / 2] = 1 - via;
	else
		(*swapbits)[(int)index / 2] = via;
	return 1 - via;
}

int SetSwapperI(vector<int>* swapbits, int index, int via)
{
	if (index % 2 == 1) {
		(*swapbits)[(int)index / 2] = 1 - via;
	}
	else
		(*swapbits)[(int)index / 2] = via;
	return 1 - via;
}

int SetSwapperO(vector<int>* swapbits, int index, int via)
{
	if (index % 2 == 0)
		(*swapbits)[(int)index / 2] = via;
	else
		(*swapbits)[(int)index / 2] = 1 - via;
	return 1 - via;
}

int neighbor(int x)
{
	if (x % 2)
		return x - 1;
	else
		return x + 1;
}

vector<int> WaksmanPermutationSetup(int* perm, int length)
{
	return sortingNetworkBits(computeInversePermutation(perm, length), length);
}

int* evaluateWaksmanNetwork_old(vector<int> swapbits, int* input, int length)
{
	if (length == 1)
		return input;
	int ins = (int)length / 2;
	int outs = (int)length / 2;
	int n = swapbits.size();
	if ((length % 2) == 0 )
		outs -= 1;
	for (int i = 0; i < ins; i++)
	{
		swapGate(input, 2 * i, 2 * i + 1, swapbits[i]);
	}
	int* input_up = (int*)malloc((int)length / 2 * sizeof(int));
	int* input_down = (int*)malloc((int)length / 2 * sizeof(int));
	for (int i = 0; i < length; i++)
	{
		if ((i % 2) == 0)
			input_up[i / 2] = input[i];
		else
			input_down[(i - 1) / 2] = input[i];
	}
	vector<int> swapbits_up(swapbits.begin() + ins, swapbits.begin() +  n/2 + 1 );
	vector<int> swapbits_down(swapbits.begin() +  n/2 + 1, swapbits.end() - outs  );

	input_up = evaluateWaksmanNetwork_old(swapbits_up, input_up, (int)length / 2);
	input_down = evaluateWaksmanNetwork_old(swapbits_down, input_down, (int)length / 2);
	for (int i = 0; i < length; i++)
	{
		if ((i % 2)==0)
			input[i] = input_up[i / 2];
		else
			input[i] = input_down[(i - 1) / 2];
	}
	for (int i = 0; i < outs; i++)
	{
		swapGate(input, 2 * i, 2 * i + 1, swapbits[n - outs + i]);
	}

	return input;
}

void swapGate(int* input, int pos1, int pos2, int swapbit)
{
	if (swapbit)
	{
		int temp = input[pos1];
		input[pos1] = input[pos2];
		input[pos2] = temp;
	}
}

vector<int> sortingNetworkBits(int* perm, int length)
{
	if (length == 1){
		vector<int> swapbits; 
		return swapbits; 
	}
	if (length == 2)
	{
		vector<int> swapbits(1);
		if (perm[0] == 0)
			swapbits[0] = 0;
		else
			swapbits[0] = 1;
		return swapbits;
	}
	int* inverse = computeInversePermutation(perm, length);
	int ins = (int)length / 2;
	int outs = (int)length / 2;
        if (length % 2 == 0) outs -=1;
	int totalSwaps = length - 1;
	vector<int> I(ins);
	vector<int> O(outs);
	for (int i = 0; i < outs; i++)
		O[i] = -1;

	for (int i = 0; i < ins; i++)
		I[i] = -1;

	int** pi = (int**)malloc(2 * sizeof(int*));
	pi[0] = (int*)malloc(((int)length / 2) * sizeof(int));
	pi[1] = (int*)malloc(((int)length / 2  + (length %2)) * sizeof(int));

	int outvia; 
	int invia; 
	if(length % 2 == 0 ){
		outvia = 1;
		int j = length - 1;
		int i = inverse[j];
		invia = SetSwapper(&I, i, outvia);
		totalSwaps -= 1;
		pi[outvia][(int)i / 2] = (int)j / 2;
		while (true)
		{
			i = neighbor(i);
			j = perm[i];
			if (j == length - 2)
			{
				pi[0][(int)i / 2] = (int)j / 2;
				break;
			}
			outvia = SetSwapper(&O, j, invia);
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
			j = neighbor(j);
			i = inverse[j];
			invia = SetSwapper(&I, i,outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
		}
	} else{
		int j = length - 1; 
		int i = inverse[j];
		outvia =1 ;
		while (i != length-1){
			invia = SetSwapper(&I, i, outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
			i = neighbor(i); 
			j = perm[i]; 
			outvia = SetSwapper(&O, j, invia); 
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
            j = neighbor(j);
            i = inverse[j];
		}
		pi[outvia][(int) i / 2] = (int)j / 2;
	}

	int j0 = outs - 1;
	while (totalSwaps > 0) {
		while (O[j0] != -1) {
			j0--;
		}
		int j = 2 * j0;
		int i = inverse[j];
		while (true) {
			invia = SetSwapper(&I, i, outvia);
			pi[outvia][(int)i / 2] = (int)j / 2;
			i = neighbor(i);
			j = perm[i];
			outvia = SetSwapper(&O, j, invia);
			pi[invia][(int)i / 2] = (int)j / 2;
			totalSwaps -= 2;
			j = neighbor(j);
			i = inverse[j];
			if (j == 2 * j0)
				break;
		}
	}

	vector<int> pi_up = sortingNetworkBits(pi[0], (int)length / 2);
	vector<int> pi_down = sortingNetworkBits(pi[1], (int)length / 2 + (length % 2));
	I.insert(I.end(), pi_up.begin(), pi_up.end());
	I.insert(I.end(), pi_down.begin(), pi_down.end());
	I.insert(I.end(), O.begin(), O.end());
	return I;
}


void evaluateWaksmanNetwork(vector<int> swapbits, int* input, int length)
{
	if (length == 1)
		return; // input;
	if (length == 2){
		swapGate(input, 0, 1, swapbits[0]);
		return;
	}
	int ins = (int)length / 2;
	// int outs = (int)length / 2-1;
	for (int i = 0; i < ins; i++)
	{
		swapGate(input, 2 * i, 2 * i + 1, swapbits[i]);
	}
	int* input_up = (int*)malloc((int)length / 2 * sizeof(int));
	int* input_down = (int*)malloc((int)((length / 2) +  (length %2)) * sizeof(int));
	for (int i = 0; i < ins; i++)
	{
		input_up[i] = input[2*i];
		input_down[i] = input[2*i+1];
	}
	if(length %2  == 1){
		input_down[ins] = input[length-1];
	}
	int n1 = count_swapbits(length / 2); 
	int n2 = count_swapbits(length / 2 + (length %2)); 
	// cout << "n1, n2 = " << n1 << ", " << n2 << endl; 


	vector<int> swapbits_up(swapbits.begin() + ins, swapbits.begin() + ins + n1);
	vector<int> swapbits_down(swapbits.begin() + ins + n1, swapbits.begin() + ins + n1 + n2);

	evaluateWaksmanNetwork(swapbits_up, input_up, (int)length / 2);
	evaluateWaksmanNetwork(swapbits_down, input_down, (int)length / 2 + (length %2));
	for (int i = 0; i < ins; i++)
	{
		input[2*i] = input_up[i];
		input[2*i + 1] = input_down[i];
	}
	if(length %2  == 1){
		input[length-1] = input_down[ins];
	}
	// for (int i = 0; i < length; i++)
	// {
	// 	if ((i % 2) == 0)
	// 		input[i] = input_up[i / 2];
	// 	else
	// 		input[i] = input_down[(i - 1) / 2];
	// }
	int outs = length /2 -1; 
	if (length %2 != 0){
		outs += 1; 
	}
	for (int i = 0; i < outs; i++)
	{	
		if (ins+n1+n2+i >= (int) swapbits.size()){
			throw invalid_argument("wtf");
		}
		swapGate(input, 2 * i, 2 * i + 1, swapbits[ins+n1+n2+i]);
	}
}
