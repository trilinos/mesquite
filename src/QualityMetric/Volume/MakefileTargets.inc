all_headers: $(QMVOLUMEHEADERS)
all_objects: $(QMVOLUMEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMVOLUMEDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QMVOLUMEDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QMVOLUMEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	