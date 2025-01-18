CXX = em++

# path to the include directory in the gemmi repo
GEMMI_DIR=../../gemmi
INCLUDE = -I$(GEMMI_DIR)/include

FLAGS = -Os -flto -Wall -Wextra -std=c++20 -s STRICT=1

LINK_FLAGS = \
  -s EXPORTED_RUNTIME_METHODS=getValue,UTF8ToString,HEAP8,HEAPF32 \
  -s EXPORTED_FUNCTIONS='[]' \
  -s INCOMING_MODULE_JS_API=onRuntimeInitialized \
  -s ALLOW_MEMORY_GROWTH=1

GEMMI_OBJS =

OBJS = density.o utils.o

density.js: $(GEMMI_OBJS) $(OBJS)
	$(CXX) $(FLAGS) $(LINK_FLAGS) $(GEMMI_OBJS) $(OBJS) -o $@

describe.h:
	echo -n '#define GEMMI_DESCRIBE "' >$@
	git -C $(GEMMI_DIR) describe | sed 's/$$/"/' >> $@

$(GEMMI_OBJS): %.o: $(GEMMI_DIR)/src/%.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $< -o $@

%.o: %.cpp describe.h
	$(CXX) $(INCLUDE) $(FLAGS) -c $<

clean:
	rm -f *.o *.js *.wasm describe.h

.PHONY: clean
