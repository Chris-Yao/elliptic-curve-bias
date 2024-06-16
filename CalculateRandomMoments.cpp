/*
CalculateRandomMoments.cpp

Calculates the moments of an elliptic curve with form y^2 = x^3 + a(t)*x^2 + b(t)*x + c(t)
for random values of t. 

Not really used during the paper; mostly as a proof of concept.

Authors: Chris Yao
*/

#include<map>
#include <iostream>
#include <cmath>
#include<vector>
#include <fstream>
#include <string>
#include <random>

using namespace std;

const __int128 MAXTRIALS = 1009; 
const __int128 MAXPRIME = 10007;
const __int128 MOMENTS = 6;
const string FAMILY = "t+2,t,0";
const string DATA_PATH = "data";


__int128 integerPower(__int128 base, __int128 exponent) {
    if (exponent == 0) 
        return 1;

    __int128 result = integerPower(base, exponent / 2);
    result *= result;

    if (exponent & 1)
        result *= base;

    return result;
}

__int128 a(__int128 t) {
    return t+2; 
}

__int128 b(__int128 t) {
    return  t;
}

__int128 c(__int128 t) {
    return 0; 
}

// printing from https://stackoverflow.com/questions/25114597/how-to-print-int128-in-g
std::ostream&
operator<<( std::ostream& dest, __int128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}



__int128 getStartingIndex(int n, vector<__int128>& primeList) {
    __int128 sum = 0;
    for (int i = 0; i < n; i++) {
        sum += primeList[i];
    }
    return sum;
}


// Legendre Symbol calcs from https://math.stackexchange.com/questions/447468/fast-legendre-symbol-calculation

int calculateMoment(__int128 maxSize, __int128 maxPower) {
    /*
    Calculates random moments of elliptic curve families

    Inputs: 
        maxSize (int): prime number to calculate the moments up to
        maxPower (int): the number of moments to calculate

    Outputs:
        A .txt file of the moments.
    */
    
    //Initializing the list of primes from .txt file
    vector<__int128> primeList;
    ifstream file(DATA_PATH + "/primes.txt");
    if (!file) {
        cout << "Error reading file!";
        return 0;
    }
    int primes;
    while (file >> primes)
    {
        if (primes > maxSize) {
            break;
        }
        primeList.push_back(primes);
    }

    // Creating array of Legendre Symbols
    int nPrimes = size(primeList);

    __int128 arr_size = getStartingIndex(nPrimes, primeList);
    __int128* legendreSymbol = new __int128[arr_size];
    for (int i = 0; i < arr_size; i++) {
        legendreSymbol[i] = -1;
    }
    int p;
    __int128 startingIndex = 0;
    for (int i = 0; i < nPrimes; i++) {
        legendreSymbol[startingIndex] = 0;
        p = primeList[i];
        int s = 0;
        for (int j = 1; j <= (p-1)/2; j++) {
            s = s + 2 * j - 1;
            if (s >= p) {
                s = s - p;
            }
            legendreSymbol[startingIndex + s] = 1;
        }
        startingIndex += primeList[i];
    }

    __int128* a_p = new __int128[arr_size];
    __int128* sumOfSymbols = new __int128[nPrimes];
    __int128 inside;
    __int128 numerator;
    __int128 correction;
    for (__int128 T = 0; T < maxSize; T++) {
        for (int j = 0; j < nPrimes; j++) {
            sumOfSymbols[j] = 0;
        }

        for (__int128 x = 0; x < maxSize; x++) {
            numerator = integerPower(x, 3) + a(T) * integerPower(x, 2) + b(T) * x + c(T);
            correction = arr_size;
            for (int i = nPrimes  - 1; primeList[i] > T && i > -1; i--) {
                p = primeList[i];
                correction -= p;
                if (p > x) {
                    inside = ((numerator % p) + p) % p; 
                    sumOfSymbols[i] += legendreSymbol[correction + inside]; 
                }
                else {
                    break;
                }
            }      
        }
        
        correction = arr_size;
        for (int i = nPrimes  - 1; primeList[i] > T && i > -1; i--) {
            p = primeList[i];
            correction -= p;
            a_p[T + correction] = sumOfSymbols[i];
        }
        if (T % (maxSize / 10 + 1) == 0) {
            cout << T / (maxSize / 10 + 1) << "\n";
        }
    }
    

    // // Writing raw data to txt in case we want it later
    // ofstream myfile (DATA_PATH + "/" + FAMILY + ",raw.txt");
    // if (myfile.is_open())
    // {
    //     for (int i = 0; i < arr_size; i++) {
    //         myfile << a_p[i] << '\n';
    //     }
    //     myfile.close();
    // }
    // else cout << "Unable to open file";

    // Writing good data to a file
    __int128 total;
    ofstream myfiletwo (DATA_PATH + "/" + FAMILY + ",random.txt");
    if (myfiletwo.is_open())
    {
        startingIndex = 0;
        for (auto &p : primeList) {
            myfiletwo << p << ',';
            for (int r = 1; r <= maxPower; r++) {
                total = 0;
                for (int count = 0; count < p; count++) {
                    total += integerPower(a_p[startingIndex + count], r);
                }
                if (r != maxPower) {
                    myfiletwo << total << ',';
                }
                else {
                    myfiletwo << total << '\n';
                }
                
            }
            startingIndex += p;
        }
        myfiletwo.close();
    }
    else cout << "Unable to open file";

    delete [] legendreSymbol;
    delete [] a_p;
    delete [] sumOfSymbols;
    return 1;
} 

int calculateRandomMoment(__int128 p) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distr(0, p-1);
    int* legendreSymbol = new int[p];
    legendreSymbol[0] = 0;
    for (int i = 1; i < p; i++) {
        legendreSymbol[i] = -1;
    }

    int s = 0;
    for (int i = 1; i <= (p-1)/2; i++) {
        s = s + 2 * i - 1;
        if (s >= p) {
            s = s - p;
        }
        legendreSymbol[s] = 1;
    }

    __int128* a_p = new __int128[MAXTRIALS];
    __int128 sumOfSymbols;
    __int128 val;
    __int128 tmp;
    __int128 T;
    for (int i = 0; i < MAXTRIALS; i++) {
        sumOfSymbols = 0;
        T = distr(gen);
        for (int x = 0; x < p; x++) {
            tmp = integerPower(x, 3) + a(T) * integerPower(x, 2) + b(T) * x + c(T);
            val = ((tmp % p) + p) % p; 
            sumOfSymbols += legendreSymbol[val];  
        }
        a_p[i] = (int) round(sumOfSymbols * (double)p / (double) MAXTRIALS);
    }

    __int128 total;
    ofstream myfiletwo (DATA_PATH + "/" + FAMILY + ",random.txt", ios_base::app);
    if (myfiletwo.is_open())
    {
        myfiletwo << p << ',';
        for (int r = 1; r <= MOMENTS; r++) {
            total = 0;
            for (int count = 0; count < MAXTRIALS; count++) {
                total += integerPower(a_p[count], r);
            }
            if (r != MOMENTS) {
                myfiletwo << total << ',';
            }
            else {
                myfiletwo << total << '\n';
            }
            
        }
        myfiletwo.close();
    }
    else {
        cout << "Unable to open file";
    }
    delete [] legendreSymbol;
    delete [] a_p;
    return 1;
}

int main(void) {
    calculateMoment(MAXTRIALS, MOMENTS);

    vector<__int128> primeList;
    ifstream file("primes.txt");
    if (!file) {
        cout << "Error reading file!";
        return 0;
    }
    int primes;
    while (file >> primes)
    {
        if (primes > MAXPRIME) {
            break;
        }
        primeList.push_back(primes);
    }

    for (int i = 0; i < size(primeList); i++) {
        if (i % (MAXTRIALS / 100 + 1) == 0) {
            cout << i / (MAXTRIALS / 100 + 1) << "\n";
        }
        if (primeList[i] < MAXTRIALS) {
            continue;
        }
        else {
            calculateRandomMoment(primeList[i]);
        }
    }
    return 0;
}