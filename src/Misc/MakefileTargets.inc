all_headers: $(MISCHEADERS)
all_objects: $(MISCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(MISCDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(MISCDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(MISCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	