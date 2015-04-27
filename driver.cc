#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector> 
#include "Hashtable.h"
#include "Unit.h" 

using namespace std;


void usage(char* argv[]) {
	cerr << "Usage: " << argv[0] << " scheme numberOfKeys sizeOfTable" << endl;
	exit(-1);
}
void unitTests( void ); 

int main(int argc, char *argv[]) {
	int n, m;
	srand(time(0));

    // run it wihout any command line parameters to invoke Unit testing
    if (argc == 1) {  
        unitTests(); 
        exit(0);  // NOTE that when unit testing,  it's all done HERE  
    }
    // or with all 4 to focus on a particular scheme 
	if (argc != 4) {
		usage(argv);
	}
	n = atoi(argv[2]);
	m = atoi(argv[3]);
	Hashtable H(m);
	
	switch(argv[1][0]) {
	case 'l':
	case 'L':
		for (int i=0; i<n; ++i) {
			H.linsert(rand() + 1);
		}
		//    H.print();
		cout << "Linear: ";  
		H.printStats();
		break;
	case 'q':
	case 'Q':
		for (int i=0; i<n; ++i) {
			H.qinsert(rand() + 1);
		}
		//    H.print();
		cout << "Quadratic: ";
		H.printStats();
		break;
	case 'd':
	case 'D':
		for (int i=0; i<n; ++i) {
			H.dinsert(rand() + 1);
		}
		//    H.print();
		cout << "Double Hashing: "; 
		H.printStats();
		break;
	default:
		usage(argv);
	}
}

// ************ DO NOT CHANGE BELOW HERE *****************
void unitTests(void) {
    Unit unit; 
    const int unitSize = 11; 
    Hashtable HT(unitSize); 
    vector<int> HV; // hashVector 
    vector<int> EV; // emptyVector
    vector<int> FV; // fullVector 
    for (int ii = 0; ii < unitSize; ++ ii) { 
       EV.push_back( EMPTY ); 
    }
    HT.getHashVector( HV ); 
    unit.assertEquals("New Hashtable _size", unitSize, HT.getSize() ); 
    unit.assertEquals("HashVector size", unitSize, HV.size() ); 
    unit.assertvectorEquals("Empty Hashtable", EV, HV );
    unit.assertEquals("New Hashtable inserts", 0, HT.getInserts() );
    unit.assertEquals("New Hashtable failures", 0, HT.getFailures() );
    unit.assertEquals("New Hashtable probes", 0, HT.getProbes() );
    // linear, no collisions 
    for (int ii = 0; ii < unitSize; ++ ii ){ 
        int key = ii + (unitSize*ii);
        HT.linsert( key ); 
        HT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "Linear with no collisions, table[" << ii << "] "; 
        unit.assertEquals(cmt.str(), key, HV[ii] ); 
        unit.assertEquals("Linear with no collisions, inserts ", ii + 1, HT.getInserts() );
        unit.assertEquals("Linear with no collisions, failures ", 0, HT.getFailures() );
        unit.assertEquals("Linear with no collisions, probes ", ii + 1, HT.getProbes() );
    }
    // linear, all collisions (table is full) 
    HT.getHashVector( FV ); 
    unit.assertvectorEquals("Linear, full table, vectors", FV, HV );
    for (int ii = 0; ii < unitSize; ++ ii ){ 
        int key = ii + (unitSize*ii);
        HT.linsert( key ); 
        HT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "Linear, all collisions, vector after insert " << ii+1 << " "; 
        unit.assertvectorEquals(cmt.str(), FV, HV );
        unit.assertEquals("Linear, all collisions, inserts ", unitSize, HT.getInserts() );
        unit.assertEquals("Linear, all collisions, failures ", ii + 1, HT.getFailures() );
        unit.assertEquals("Linear, all collisions, probes ", unitSize, HT.getProbes() );
    }    
    // clear table for another test 
    HT.clear(); 
    HT.getHashVector( HV ); 
    unit.assertEquals("Cleared Hashtable _size", unitSize, HT.getSize() ); 
    unit.assertvectorEquals("Cleared Hashtable", EV, HV );
    unit.assertEquals("Cleared Hashtable inserts", 0, HT.getInserts() );
    unit.assertEquals("Cleared Hashtable failures", 0, HT.getFailures() );
    unit.assertEquals("Cleared Hashtable probes", 0, HT.getProbes() );
    // linear, some collisions 
    HT.getHashVector( FV ); // put expected values in FV  
    HT.linsert( unitSize * unitSize ); 
    FV[0] = unitSize * unitSize; 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 0", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 0, inserts=", 1, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 0, probes=", 1, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 0, failures=", 0, HT.getFailures() );
    HT.linsert( 1 + unitSize );
    FV[1] =  1 + unitSize; 
    HT.linsert( 1 + (unitSize*unitSize) );
    FV[2] =  1 + (unitSize*unitSize); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 1,2", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 1,2 inserts=", 3, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 1,2 probes=", 4, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 1,2 failures=", 0, HT.getFailures() );
    HT.linsert( 2 + unitSize );
    FV[3] =  2 + unitSize; 
    HT.linsert( 2 + (unitSize*unitSize) );
    FV[4] =  2 + (unitSize*unitSize); 
    HT.linsert( 2 + unitSize );
    FV[5] =  2 + unitSize; 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 3,4,5", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 3,4,5 inserts=", 6, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 3,4,5 probes=", 13, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 3,4,5 failures=", 0, HT.getFailures() );
    HT.linsert( 6 );
    FV[6] =  6; 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 6", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 6 inserts=", 7, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 6 probes=", 14, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 6 failures=", 0, HT.getFailures() );
    HT.linsert( 10 + unitSize );
    FV[10] = 10  + unitSize; 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 10", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 10 inserts=", 8, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 10 probes=", 15, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 10 failures=", 0, HT.getFailures() );
    HT.linsert( 10 + (unitSize*unitSize) );
    FV[7] = 10  + (unitSize*unitSize); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 7", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 7 inserts=", 9, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 7 probes=", 24, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 7 failures=", 0, HT.getFailures() );
    HT.linsert( 10 + (unitSize*unitSize*unitSize) );
    FV[8] = 10  + (unitSize*unitSize*unitSize); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 8", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 8 inserts=", 10, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 8 probes=", 34, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 8 failures=", 0, HT.getFailures() );
    HT.linsert( 9 );
    FV[9] = 9;
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, insert at 9", FV, HV );
    unit.assertEquals("Linear, some collisions, insert at 8 inserts=", 11, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, insert at 8 probes=", 35, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, insert at 8 failures=", 0, HT.getFailures() );
    HT.linsert( 8 + unitSize );
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Linear, some collisions, failure", FV, HV );
    unit.assertEquals("Linear, some collisions, failure, inserts=", 11, HT.getInserts() );
    unit.assertEquals("Linear, some collisions, failure, probes=", 35, HT.getProbes() );
    unit.assertEquals("Linear, some collisions, failure, failures=", 1, HT.getFailures() );
    
    
    // clear table for another test 
    HT.clear(); 
    // quadratic, no collisions 
    for (int ii = 0; ii < unitSize; ++ ii ){ 
        int key = ii + (unitSize*ii);
        HT.qinsert( key ); 
        HT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "Quadratic with no collisions, table[" << ii << "] "; 
        unit.assertEquals(cmt.str(), key, HV[ii] ); 
        unit.assertEquals("Quadratic with no collisions, inserts ", ii + 1, HT.getInserts() );
        unit.assertEquals("Quadratic with no collisions, failures ", 0, HT.getFailures() );
        unit.assertEquals("Quadratic with no collisions, probes ", ii + 1, HT.getProbes() );
    }
    // quadratic, all collisions (table is full) 
    HT.getHashVector( FV ); 
    unit.assertvectorEquals("Quadratic, full table, vectors", FV, HV );
    for (int ii = 0; ii < unitSize; ++ ii ){ 
        int key = ii + (unitSize*ii);
        HT.qinsert( key ); 
        HT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "Quadratic, all collisions, vector after insert " << ii+1 << " "; 
        unit.assertvectorEquals(cmt.str(), FV, HV );
        unit.assertEquals("Quadratic, all collisions, inserts ", unitSize, HT.getInserts() );
        unit.assertEquals("Quadratic, all collisions, failures ", ii + 1, HT.getFailures() );
        unit.assertEquals("Quadratic, all collisions, probes ", unitSize, HT.getProbes() );
    }   
    // clear for another test  
    HT.clear(); 
    HT.getHashVector( FV ); 
    HT.qinsert( 0 + (unitSize * 1)); 
    FV[0] =   ( 0 + (unitSize * 1)); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 1st, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 1st, inserts", 1, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 1st, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 1st, probes", 1, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 2)); 
    FV[1] =   ( 0 + (unitSize * 2)); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 2nd, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 2nd, inserts", 2, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 2nd, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 2nd, probes", 3, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 3)); 
    FV[4] =   ( 0 + (unitSize * 3)); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 3rd, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 3rd, inserts", 3, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 3rd, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 3rd, probes", 6, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 4)); 
    FV[9] =   ( 0 + (unitSize * 4)); 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 4th, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 4th, inserts", 4, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 4th, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 4th, probes", 10, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 5)); 
    FV[5] =   ( 0 + (unitSize * 5)); // 16 % 11 = 5
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 5th, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 5th, inserts", 5, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 5th, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 5th, probes", 15, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 6)); 
    FV[3] =   ( 0 + (unitSize * 6)); // 25 % 11 = 3
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 6th, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 6th, inserts", 6, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 6th, failures", 0, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 6th, probes", 21, HT.getProbes() );
    HT.qinsert( 0 + (unitSize * 7)); 
//  FV[??] =   ( 0 + (unitSize * 7)); // 36 % 11 = 3; 49 % 11 = 5; 64 % 11 = 9; 81 % 11 = 4; 100 % 11 = 1 
    HT.getHashVector( HV ); 
    unit.assertvectorEquals("Quadratic, some collisions, after 7th, vectors ", FV, HV );
    unit.assertEquals("Quadratic, some collisions, after 7th, inserts", 6, HT.getInserts() );
    unit.assertEquals("Quadratic, some collisions, after 7th, failures", 1, HT.getFailures() );
    unit.assertEquals("Quadratic, some collisions, after 7th, probes", 21, HT.getProbes() );

    // clear table for another test 
    const int dhSize = 37;
    Hashtable DT(37);
    // double, no collisions 
    for (int ii = 0; ii < dhSize; ++ ii ){ 
        int key = ii + (dhSize*ii);
        DT.dinsert( key ); 
        DT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "DoubleHashing with no collisions, table[" << ii << "] "; 
        unit.assertEquals(cmt.str(), key, HV[ii] ); 
        unit.assertEquals("DoubleHashing with no collisions, inserts ", ii + 1, DT.getInserts() );
        unit.assertEquals("DoubleHashing with no collisions, failures ", 0, DT.getFailures() );
        unit.assertEquals("DoubleHashing with no collisions, probes ", ii + 1, DT.getProbes() );
    }
    // quadratic, all collisions (table is full) 
    DT.getHashVector( FV ); 
    unit.assertvectorEquals("DoubleHashing, full table, vectors", FV, HV );
    for (int ii = 0; ii < dhSize; ++ ii ){ 
        int key = ii + (dhSize*ii);
        DT.dinsert( key ); 
        DT.getHashVector( HV ); 
        std::stringstream cmt; 
        cmt << "DoubleHashing, all collisions, vector after insert " << ii+1 << " "; 
        unit.assertvectorEquals(cmt.str(), FV, HV );
        unit.assertEquals("DoubleHashing, all collisions, inserts ", dhSize, DT.getInserts() );
        unit.assertEquals("DoubleHashing, all collisions, failures ", ii + 1, DT.getFailures() );
        unit.assertEquals("DoubleHashing, all collisions, probes ", dhSize, DT.getProbes() );
    }   
    
    unit.printResults();


} 


