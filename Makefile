hash: Hashtable.cc Hashtable.h Unit.cc Unit.h driver.cc
	g++ -Wall -g -c Unit.cc -o Unit.o 
	g++ -Wall -g -c Hashtable.cc -o Hashtable.o 
	g++ -Wall -g driver.cc Hashtable.o Unit.o -o hash

clean: 
	rm -f hash
	rm -f HashTable.o 
	rm -f Unit.o 
	rm -fr hash.dSYM/* 
	rmdir hash.dSYM
 
