#include <cstdlib>
#include <iostream>
#include <cassert> 
#include <vector> 
#include "Hashtable.h"

using namespace std; 

#define HUGE_FLOAT 9999999.9999999 

Hashtable::Hashtable(int size) {
	//constructor
	int prime = nextPrime( size ); 
	_size = size;
	if (prime != size) { 
	  cout << "Warning: size = " << size << " is not a prime number."; 
//	  cout << "Using " << prime << " instead." << endl; 
//	  _size = prime; 
	}
	_totalProbes = 0;
	_numInserts = 0;
	_numFailures = 0; 
	_table = new int [ _size ];
	for (int ii=0; ii < _size; ii++) {
		_table[ ii ] = EMPTY;
	}
}
void Hashtable::clear(void) {
	_totalProbes = 0;
	_numInserts = 0;
	_numFailures = 0; 
	for (int ii=0; ii < _size; ii++) {
		_table[ ii ] = EMPTY;
	}
}    
    

Hashtable::~Hashtable() {
	//deconstructor
	delete[] _table;
}

void Hashtable::tallyProbes(int p) {
	// Add the number of probes, p, for one insert operation to the tally.
	_totalProbes += p;
	_numInserts += 1;
}

void Hashtable::printStats() {
	cout << "Average probes/insert = " <<
			probeRate() << " = " <<
			_totalProbes << "/" << _numInserts;
	cout <<", size = "<< _size;
	cout <<", failures = " << _numFailures << endl;
}

float Hashtable::probeRate() {
    if ( _numInserts == 0 ) return HUGE_FLOAT;  
	return (float)_totalProbes / (float)_numInserts;
}

int Hashtable::hash(int k) {
	return k % _size;
}

int Hashtable::hash2(int k, int iter) {
	// ******* CHANGE THIS ***************
	return 31 - (k % 31);
	
}

void Hashtable::qinsert(int k) {
	// Insert k in the hash table.
	// Use open addressing with quadratic probing and hash(k) = k % _size.

	// Tips:
	// - Look at the lecture slides if you don't remember how to implement these.
	// - You need to prevent against infinite loops. You should limit the number
	// of times you probe and print an error message upon exceeding this limit.
	// - You're welcome to add new declarations to hash.h, add to the constructor,
	// include libraries, or anything else you need. Just make it work.
	// - You're also welcome to modify the main() method to automate your testing.
	// In fact, you should!

    // ************* ADD YOUR CODE HERE *******************
	if (_numInserts == _size) {
		cout << "Warning: qinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl;
		_numFailures++;
		return;
	}
	else{
		int probepoint_q = hash(k);
	    int i = 1;
		int counts_q=1;
		if (_table[probepoint_q]==EMPTY) {   //no tallyProbes() because there is no probing 
			_table[probepoint_q]=k;
			tallyProbes(counts_q);
			return;
		}
		else {
			for (int q = 0 ; q < _size; q++){
                probepoint_q=(probepoint_q+2*i-1)%_size;
				counts_q++;
				if (_table[probepoint_q]==EMPTY) {
					_table[probepoint_q]=k;
					tallyProbes(counts_q);
					return;
				}
		        i++;
			}
		}
		//tallyProbes(counts_q+1);
		cout << "Warning: qinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl;
		_numFailures++;
		return;
	}
    // Your method should return after it stores the value in an EMPTY slot 
    // and calls tallyProbes, so if it gets here, it didn't find an EMPTY slot 
    //_numFailures += 1; 
    cout << "Warning: qinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl;

    	
}

void Hashtable::linsert(int k) {
	// Insert k in the hash table.
	// Use open addressing with linear probing and hash(k) = k % _size.

    // ************* ADD YOUR CODE HERE *******************
	if (_numInserts == _size) {
		 cout << "Warning: linsert(" << k << ") found no EMPTY slot, so no insert was done." << endl; 
		 _numFailures++;
		return;
		}
	else{
		int counts = 0;
		int probepoint_l = hash(k);
	    while (_table[probepoint_l]!=EMPTY){
		    probepoint_l=(probepoint_l+1)%_size;
			counts++;
	    }
	    _table[probepoint_l]=k;
	    tallyProbes(1+counts);
	    return;
	}
    // Your method should return after it stores the value in an EMPTY slot 
    // and calls tallyProbes, so if it gets here, it didn't find an EMPTY slot 
    //_numFailures += 1; 
 
}

void Hashtable::dinsert(int k) {
	// Insert k in the hash table.
	// Use open addressing with double hashing. 
	// Use the existing hash(int) function
	// and implement hash2(int, int) as well.

    // ************* ADD YOUR CODE HERE *******************
	if (_numInserts == _size){
		cout << "Warning: dinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl;
		_numFailures++;
		return;
	} 
	else{
		int probepoint_f= hash(k);
		int probepoint_d=0;
		int counts_d=1; //numbers of failure probing 
	    if (_table[probepoint_f]==EMPTY) { //no tallyProbes() because there is no probing 
			_table[probepoint_f]=k;
			tallyProbes(counts_d);
			return;
		}
		else{
			for (int iter = 0; iter < _size; iter++){
				 probepoint_d=(probepoint_f+iter*hash2(k,iter)) % _size;
				 counts_d++;
				 if (_table[probepoint_d]==EMPTY){
				    _table[probepoint_d]=k;
				    tallyProbes(counts_d+1);
				    return;
				 }
			}
		tallyProbes(counts_d+1);
		cout << "Warning: dinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl;
		_numFailures++;
		return;
		}
}
    // Your method should return after it stores the value in an EMPTY slot 
    // and calls tallyProbes, so if it gets here, it didn't find an EMPTY slot 
    //_numFailures += 1; 
    cout << "Warning: dinsert(" << k << ") found no EMPTY slot, so no insert was done." << endl; 
}

void Hashtable::print() {
	// Print the content of the hash table.

	for (int i=0; i<_size; i++) {
		cout<<"["<<i<<"] ";
		if (_table[i] != EMPTY)
			cout << _table[i];
		cout << endl;
	}
}

// 

int Hashtable::nextPrime( int n ) { 
    int loops = (n < 100) ? 100 : n;
    for (int ii = 0; ii < loops; ii++ ){ 
      if ( isPrime( n + ii ) ) return (n + ii); 
    }
    assert( false ); // logic error 
}
bool Hashtable::isPrime( int n ) { 
    if (n < 2) return false; 
    if ((n == 2) || (n == 3) || (n == 5)) return true;
    if (n % 2 == 0) return false;  
    if (n % 3 == 0) return false;  
    if (n % 5 == 0) return false;  
    return isPrime( n, 7 ); 
}
bool Hashtable::isPrime( int n, int divisor ){ 
    if ((divisor * divisor) > n) return true; 
    if ((n % divisor) == 0) return false; 
    assert( (divisor % 2) == 1 ); // divisor is odd  
    if ((divisor % 10) == 3) return isPrime( n, divisor + 4 ); 
    return isPrime( n, divisor + 2 ); 
}
void Hashtable::getHashVector(std::vector<int>& hashVector ){ 
    hashVector.clear(); 
    for (int ii=0; ii < _size; ii++) { 
        hashVector.push_back( _table[ii] ); 
    }
} 
int Hashtable::getSize() { 
    return _size; 
} 
int Hashtable::getInserts(){
    return _numInserts; 
}
int Hashtable::getFailures(){
    return _numFailures; 
} 
int Hashtable::getProbes(){
    return _totalProbes;
}


