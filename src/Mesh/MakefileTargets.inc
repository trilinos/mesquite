all_headers: $(MESHHEADERS)
all_objects: $(MESHOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(MESHDIR)/%.hpp
	${PREFIX} ln -sf ../$< $@

$(localincludedir)/%.hh : $(MESHDIR)/%.hh
	${PREFIX} ln -sf ../$< $@

$(localobjdir)/%.o : $(MESHDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

