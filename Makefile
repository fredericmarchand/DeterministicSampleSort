all:	createDirs main

main:   
		mpiCC main.cpp -o main

createDirs:
			mpiCC createTempDirs.cpp -o createDirs

clean:
		rm -rf main createDirs
