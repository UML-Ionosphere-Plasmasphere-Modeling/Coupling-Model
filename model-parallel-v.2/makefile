objects= main.o fieldsgrids.o particles.o parameters.o module.o module_0.o module_1.o
program: $(objects)
	h5c++ -o program $(objects) -fopenmp -O0
main.o: main.cpp particles.h parameters.h vector3.h module.h
	h5c++ -c -g -fopenmp main.cpp -O0
parameters.o: parameters.cpp parameters.h
	h5c++ -c -g -fopenmp parameters.cpp -O0
module.o: module.cpp parameters.h fieldsgrids.h particles.h vector3.h module.h
	h5c++ -c -g -fopenmp module.cpp -O0
module_0.o: module_0.cpp parameters.h fieldsgrids.h particles.h vector3.h module_0.h
	h5c++ -c -g -fopenmp module_0.cpp -O0
module_1.o: module_1.cpp parameters.h fieldsgrids.h particles.h vector3.h module_1.h
	h5c++ -c -g -fopenmp module_1.cpp -O0
fieldsgrids.o: fieldsgrids.cpp fieldsgrids.h vector3.h parameters.h
	h5c++ -c -g -fopenmp fieldsgrids.cpp -O0
particles.o: particles.cpp particles.h mathutil.h vector3.h parameters.h structdef.h
	h5c++ -c -g -fopenmp particles.cpp -O0
.PHONY: clean
clean:
	rm program $(objects)
