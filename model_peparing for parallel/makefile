objects= main.o fieldsgrids.o particles.o parameters.o module.o module_0.o module_1.o
program: $(objects)
	h5c++ -o program $(objects) -fopenmp 
main.o: main.cpp particles.h parameters.h vector3.h module.h
	h5c++ -c -g -fopenmp main.cpp
parameters.o: parameters.cpp parameters.h
	h5c++ -c -g -fopenmp parameters.cpp
module.o: module.cpp parameters.h fieldsgrids.h particles.h vector3.h module.h
	h5c++ -c -g -fopenmp module.cpp
module_0.o: module_0.cpp parameters.h fieldsgrids.h particles.h vector3.h module_0.h
	h5c++ -c -g -fopenmp module_0.cpp
module_1.o: module_1.cpp parameters.h fieldsgrids.h particles.h vector3.h module_1.h
	h5c++ -c -g -fopenmp module_1.cpp
fieldsgrids.o: fieldsgrids.cpp fieldsgrids.h vector3.h parameters.h
	h5c++ -c -g -fopenmp fieldsgrids.cpp
particles.o: particles.cpp particles.h mathutil.h vector3.h parameters.h structdef.h
	h5c++ -c -g -fopenmp particles.cpp
.PHONY: clean
clean:
	rm program $(objects)
