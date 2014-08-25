icpc -c ../lib/*.c
icpc -O2 -xAVX ../gmxstress.c *.o -o stress_intel.exe -lxdrfile -L../build
rm *.o
