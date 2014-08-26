# Compile script for the GNU version of the stress code.
g++ -c ../lib/*.c
g++ -O ../gmxstress.c *.o -o stress_gnu.exe -lxdrfile -L../build
rm *.o
