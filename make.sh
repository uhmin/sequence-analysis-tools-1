mkdir -p bin
echo "compiling orf.exe"
gcc -Wall -O3 -D_FILE_OFFSET_BITS=64 -o bin/orf.exe orf.091208.c

echo "compiling translate.exe"
gcc -Wall -O3 -D_FILE_OFFSET_BITS=64 -o bin/translate.exe translate.091221.c

echo "compiling revcmp.exe"
gcc -Wall -O3 -D_FILE_OFFSET_BITS=64 -o bin/revcmp.exe revcmp2.c

echo "compiling docono.exe"
gcc -Wall -O3 -D_FILE_OFFSET_BITS=64 -o bin/docono.exe docono3.c

echo "compiling filterN.exe"
#gcc -Wall -lm -O -D_FILE_OFFSET_BITS=64 -o bin/filterN.exe filterN_090601.c
gcc -Wall -lm -O -D_FILE_OFFSET_BITS=64 -o bin/filterN.exe filterN_140630.c

echo "all compile done!"
