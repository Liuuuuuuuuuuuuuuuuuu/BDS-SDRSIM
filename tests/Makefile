CC=gcc
CFLAGS?=-I.. -Iunity -Iext/kissfft -O2 -Wall

TESTS= test_prn.o test_crc.o test_navframe.o test_chipseq.o        test_scaling.o test_fft.o helpers.o unity/unity.o        ext/kissfft/kiss_fft.o
PROJ= ../globals.o ../bch.o ../navbits.o ../channel.o ../rinex.o       ../orbits.o ../coord.o ../path.o ../bdssim.o

all: run_all

run_all: $(TESTS) run_all.o $(PROJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

../%.o: ../%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o run_all ../*.o

coverage:
	$(MAKE) CFLAGS+="--coverage -O0" all

.PHONY: all clean coverage
