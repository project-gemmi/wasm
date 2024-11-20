CXX = em++

# path to the include directory in the gemmi repo
GEMMI_DIR=../../gemmi
INCLUDE = -I$(GEMMI_DIR)/include

FLAGS = \
  -Os -flto -Wall -Wextra -std=c++20 \
  -s STRICT=1 \
  -s DISABLE_EXCEPTION_CATCHING=0

LINK_FLAGS = \
  -s EXPORTED_RUNTIME_METHODS=ccall,writeArrayToMemory,getValue,UTF8ToString,HEAP8 \
  -s EXPORTED_FUNCTIONS=_malloc \
  -s INCOMING_MODULE_JS_API=print,printErr,setStatus,onRuntimeInitialized \
  -s ALLOW_MEMORY_GROWTH=1

GEMMI_OBJS = eig3.o intensit.o mmcif.o mtz.o mtz2cif.o polyheur.o pdb.o \
	     resinfo.o to_mmcif.o to_pdb.o sprintf.o xds_ascii.o

OBJS = to_cif.o from_cif.o utils.o

convert.js: $(GEMMI_OBJS) $(OBJS)
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
