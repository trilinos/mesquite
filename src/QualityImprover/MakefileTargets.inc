all_headers: $(QUALIMPROVERHEADERS)
all_objects: $(QUALIMPROVEROBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALIMPROVERDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QUALIMPROVERDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QUALIMPROVERDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	