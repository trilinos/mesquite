all_headers: $(VERTEXMOVERHEADERS)
all_objects: $(VERTEXMOVEROBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(VERTEXMOVERDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(VERTEXMOVERDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(VERTEXMOVERDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	