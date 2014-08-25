g++ -c ../lib/*.c
#libtool --mode=link g++ ../gmxstress.c *.o -o stress_gnu.exe libxdrfile.la
g++ -O ../gmxstress.c *.o -o stress_gnu.exe -lxdrfile -L../build
rm *.o
#-L /home/afenley/Research/Gilson/Projects/ccode/libxdrfile
