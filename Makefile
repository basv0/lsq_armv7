
CC=g++
FLAGS=-O2
FLAGSNEON=-O2 -mfpu=neon
LIB=-lrt

all: mnk_test_optim mnk_test_neonh mnk_test_neonv

mnk_test_optim: mnk_test_optim.cpp
	$(CC) -o mnk_test_optim mnk_test_optim.cpp $(FLAGS) $(LIB)

mnk_test_neonh: mnk_test_neonh.cpp
	$(CC) -o mnk_test_neonh mnk_test_neonh.cpp $(FLAGSNEON) $(LIB)
	
mnk_test_neonv: mnk_test_neonv.cpp
	$(CC) -o mnk_test_neonv mnk_test_neonv.cpp $(FLAGSNEON) $(LIB)
