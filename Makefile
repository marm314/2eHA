# Makefile for 2eHA
#
#CPP = icc
#CPPFLAGS = -Wall 
CPP = g++
CPPFLAGS = -Wall -O2
Cln = /bin/rm -rf
NAME=2eHA
###########################################
###########################################
SCR= main.cpp D_read_calc_rho.cpp HA_sol_n.cpp Input_commands.cpp legendre_quadrature.cpp Mathematical_Functions.cpp sphere_lebedev_rule.cpp String_ops.cpp 
OBJECTS= main.o D_read_calc_rho.o HA_sol_n.o Input_commands.o legendre_quadrature.o Mathematical_Functions.o sphere_lebedev_rule.o String_ops.o
LIB=libcuba.a

all:
	make 2ehasol

2ehasol: $(OBJECTS) $(SCR)  $(LIB) cuba.h Makefile 
	$(CPP) $(CPPFLAGS) $(OBJECTS) $(LIB) -o $(NAME)_EXE

%.o: %.cpp   
	$(CPP) $(CPPFLAGS) -c $*.cpp 

clean:
	$(Cln) *.o
	$(Cln) *$(NAME)_EXE
	$(Cln) *~
	$(Cln) $(NAME).tar.gz 
tar:
	mkdir $(NAME)
	cp *.cpp *.h Makefile ./$(NAME)
	tar -pczf $(NAME).tar.gz ./$(NAME)
	rm -r ./$(NAME)
	cp ./*.tar.gz ../
