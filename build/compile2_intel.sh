# Compile script using Intel compilers (if available). Make sure the 
# library file libxdrfile.a was also compiled using the Intel compilers!
icpc -c ../lib/*.c
icpc -O2 -xAVX ../gmxstress.c *.o -o stress_intel.exe -lxdrfile -L../build
rm *.o
