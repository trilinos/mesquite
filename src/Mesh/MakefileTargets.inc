all_headers: $(MESHHEADERS)
all_objects: $(MESHOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(MESHDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(MESHDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(MESHDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	