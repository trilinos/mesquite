all_headers: $(RANDSMOOTHHEADERS)
all_objects: $(RANDSMOOTHOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(RANDSMOOTHDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(RANDSMOOTHDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(RANDSMOOTHDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

