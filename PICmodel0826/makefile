objects= main.o fieldsgrids.o particles.o parameters.o module.o
program: $(objects)
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
.PHONY: clean
clean:
	rm program $(objects)
