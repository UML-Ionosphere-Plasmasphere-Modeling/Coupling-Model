objects= main.o fieldsgrids.o particles.o parameters.o module.o
program: $(objects)
<<<<<<< HEAD:PICmodel/makefile
	h5c++ -o program $(objects)
main.o: main.cpp particles.h parameters.h vector3.h module.h
	h5c++ -c -g main.cpp
parameters.o: parameters.cpp parameters.h
	h5c++ -c -g parameters.cpp
module.o: module.cpp parameters.h fieldsgrids.h particles.h vector3.h module.h
	h5c++ -c -g module.cpp
fieldsgrids.o: fieldsgrids.cpp fieldsgrids.h vector3.h parameters.h
	h5c++ -c -g fieldsgrids.cpp
particles.o: particles.cpp particles.h mathutil.h vector3.h parameters.h structdef.h
	h5c++ -c -g particles.cpp
=======
	h5c++ -o program $(objects) -fopenmp 
main.o: main.cpp particles.h parameters.h vector3.h module.h
	h5c++ -c -g -fopenmp main.cpp
parameters.o: parameters.cpp parameters.h
	h5c++ -c -g -fopenmp parameters.cpp
module.o: module.cpp parameters.h fieldsgrids.h particles.h vector3.h module.h
	h5c++ -c -g -fopenmp module.cpp
fieldsgrids.o: fieldsgrids.cpp fieldsgrids.h vector3.h parameters.h
	h5c++ -c -g -fopenmp fieldsgrids.cpp
particles.o: particles.cpp particles.h mathutil.h vector3.h parameters.h structdef.h
	h5c++ -c -g -fopenmp particles.cpp
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/makefile
.PHONY: clean
clean:
	rm program $(objects)
