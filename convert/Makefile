CXX = em++

# path to the include directory in the gemmi repo
INCLUDE = -I../../gemmi/include

FLAGS = \
  -Os -flto -Wall -Wextra -std=c++20 \
  -s WASM=1 \
  -s STRICT=1 \
  -s DISABLE_EXCEPTION_CATCHING=0 \
  -s ALLOW_MEMORY_GROWTH=1

LINK_FLAGS = \
  -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["ccall", "writeArrayToMemory"]' \
  -s 'EXPORTED_FUNCTIONS=[]'

OBJS = to_cif.o utils.o write.o cif2pdb.o

convert.js: $(OBJS)
	$(CXX) $(FLAGS) $(LINK_FLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $<

clean:
	rm -f *.o *.js *.wasm

.PHONY: clean
