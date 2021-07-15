#Wojciech :Krzemien 13.02.2006
#Makefile for DriftChamber simulation
# 
# make run, make clean, make tar
#===============================================================================
DIR=tragaldabas
###################
ARCH         := $(shell root-config --arch)
PLATFORM      = $(ARCH)
###################
NAME1 = tunpacker
NAME2 = thldevent
NAME3 = thldsubevent
NAME4 = tevent
NAME5 = teventhdr
NAME6 = tmydict
NAME7 = thit
NAME8 = trpclookuptable
NAME9 = trpccalpar
NAME10 = trpcraw
NAME11 = trpcrawf
NAME12 = trpchit
NAME13 = trpchitf
NAME14 = trpcsaeta
NAME15 = trpcsaetaf
NAME16 = tactivecells
NAME17 = ttmatrix
###################
LINKDEF = LinkDef.h
###################
EXEC1 = myroot.x
###################
OBJS1 = $(NAME1).o
OBJS2 = $(NAME2).o
OBJS3 = $(NAME3).o
OBJS4 = $(NAME4).o
OBJS5 = $(NAME5).o
OBJS6 = $(NAME6).o
OBJS7 = $(NAME7).o
OBJS8 = $(NAME8).o
OBJS9 = $(NAME9).o
OBJS10 = $(NAME10).o
OBJS11 = $(NAME11).o
OBJS12 = $(NAME12).o
OBJS13 = $(NAME13).o
OBJS14 = $(NAME14).o
OBJS15 = $(NAME15).o
OBJS16 = $(NAME16).o
OBJS17 = $(NAME17).o
###################
HEADS1 = $(NAME1).h 
HEADS2 = $(NAME2).h 
HEADS3 = $(NAME3).h 
HEADS4 = $(NAME4).h 
HEADS5 = $(NAME5).h 
HEADS6 = $(NAME6).h 
HEADS7 = $(NAME7).h 
HEADS8 = $(NAME8).h 
HEADS9 = $(NAME9).h 
HEADS10 = $(NAME10).h 
HEADS11 = $(NAME11).h 
HEADS12 = $(NAME12).h 
HEADS13 = $(NAME13).h 
HEADS14 = $(NAME14).h 
HEADS15 = $(NAME15).h 
HEADS16 = $(NAME16).h
HEADS17 = $(NAME17).h   
###################
LIB_SHARED = lib$(NAME1).so
###################
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
###################
CFLAGS = -O -Wall -g -fPIC
LDFLAGS=-O
###################
CFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
###################
CO = g++
LD = $(CO)
###################
%.o: %.cc %.h
	$(CO) $(CFLAGS) -c $<
%.o: %.cc
	$(CO) $(CFLAGS) -c $<
###################
$(LIB_SHARED): $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS8) $(OBJS7) $(OBJS9) $(OBJS10) $(OBJS11) $(OBJS12) $(OBJS13) $(OBJS14) $(OBJS15) $(OBJS16) $(OBJS17)
	$(LD) -shared $(LDFLAGS) $^ $(LIBS) $(GLIBS) -o $@
###################
run: $(EXEC1)
	./$(EXEC1)
###################
clean:
	rm -f *.o *.x core core* a.out *.so tmydict.*
###################
tar: clean
	(cd ../; tar -cvzf $(DIR).tar.gz $(DIR) )
###################
$(NAME6).cc:$(HEADS1) $(HEADS2) $(HEADS3) $(HEADS4) $(HEADS5) $(HEADS8) $(HEADS7) $(HEADS9) $(HEADS10) $(HEADS11) $(HEADS12) $(HEADS13) $(HEADS14) $(HEADS15)  $(HEADS16) $(HEADS17)    $(LINKDEF)
	@echo "Generation dictrionary................................................."
	@rootcint -f $(NAME6).cc -c -p $(HEADS1) $(HEADS2) $(HEADS3) $(HEADS4) $(HEADS5) $(HEADS8) $(HEADS7) $(HEADS9) $(HEADS10) $(HEADS11) $(HEADS12) $(HEADS13) $(HEADS14) $(HEADS15)  $(HEADS16) $(HEADS17)  $(LINKDEF)
###################
