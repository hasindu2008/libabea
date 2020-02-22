CC       = gcc
CXX      = g++
CFLAGS   += -g  -Wall -O2 -std=c++11 
LDFLAGS  += $(LIBS) 
BUILD_DIR = .
BINARY = abea
OBJ = main.o \
	  f5c.o \
      events.o \
      model.o \
      align.o 


$(BINARY): $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: main.c f5c.h 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/f5c.o: f5c.c f5c.h 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/events.o: events.c f5c.h  ksort.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/model.o: model.c model.h f5c.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/align.o: align.c f5c.h 
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@


clean: 
	rm -rf $(BINARY) $(BUILD_DIR)/*.o


