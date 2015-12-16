
CC        = clang
LDFLAGS   = -O 
CFLAGS    = -Weverything -Wextra -pedantic $(LDFLAGS)
LDLIBS    = $(shell gsl-config --libs)

.SUFFIXES:
.SUFFIXES:  .c .o .h

.PHONY: clean indent

target    = vegas-integral

$(target): $(target).o two-cubes-integrand2.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

png: res
	gnuplot two-cubes.gp
	
res: $(target)
	./$(target) > res

indent : $(target).c
	indent $<

clean : 
	rm -f *.o
	rm -f *.*~
	rm  -f $(target)
