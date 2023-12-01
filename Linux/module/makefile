CC = g++

SOURCES = mmModules.cpp

LIBS = -lm -lpthread
PROGRAM = mmModules

MAKEDEP = domakedep

OBJECTS= $(SOURCES:.c=.o)
LINTFILES=$(SOURCES:.c=.ln)

CFLAGS = 
# -g

$(PROGRAM): $(OBJECTS) 
	$(CC) $(CFLAGS) -O3 -o $@ $(OBJECTS) $(LIBS)

$(LINTFILES):
	$(LINT) $(LINTFILES) 

clean:
	rm  compos2.o

#clall : clean all
clall:
	rm -f *.o core a.out $(PROGRAM)

depend:
	$(MAKEDEP) *.[chC] > deps_top

#include deps_top
# DO NOT DELETE THIS LINE -- ccdep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.



# IF YOU PUT ANYTHING HERE IT WILL GO AWAY

