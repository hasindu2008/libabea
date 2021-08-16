CC       = gcc
CXX      = g++
CFLAGS   += -g  -Wall -O2 -std=c++11
LDFLAGS  += $(LIBS)
BUILD_DIR = .
BINARY = abea_example

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

OBJ =   example.o \
        libabea.o \
        f5c.o \
        events.o \
        model.o \
        align.o

$(BINARY): $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILD_DIR)/example.o: example.cpp f5c.h example.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/libabea.o: libabea.cpp libabea.h f5c.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/f5c.o: f5c.cpp f5c.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/events.o: events.cpp f5c.h  ksort.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/model.o: model.cpp model.h f5c.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/align.o: align.cpp f5c.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

pylib:
	python3 setup.py build && cp build/lib.*/*.so  ./

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o
	python3 setup.py clean
	rm -rf build abea.cpp
	rm -f *.so
	rm -f out.txt outpy.txt
