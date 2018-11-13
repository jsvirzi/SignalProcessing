all: lib obj lib/libfilters.so

lib: src/filters.cpp
	mkdir -p lib

obj: src/filters.cpp
	mkdir -p obj

lib/libfilters.so: obj/filters.o
	ld -shared -soname $@ -o $@ obj/filters.o -lc -lm

obj/filters.o: src/filters.cpp include/filters.h
	gcc -fPIC -c -Iinclude $< -o obj/filters.o

examples: examples/lowpass examples/minimal
	@echo "examples complete"

examples/lowpass:
	gcc -o examples/lowpass -Iinclude examples/lowpass.cpp -lm -Llib -lfilters

examples/minimal:
	gcc -o examples/minimal -Iinclude examples/minimal.cpp -lm -Llib -lfilters

clean:
	rm -f lib/*.so obj/*.o examples/lowpass examples/minimal

