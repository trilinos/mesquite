all_headers: $(CONJGRADHEADERS)
all_objects: $(CONJGRADOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(CONJGRADDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(CONJGRADDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(CONJGRADDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	