all_headers: $(OBJFUNCHEADERS)
all_objects: $(OBJFUNCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(OBJFUNCDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(OBJFUNCDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(OBJFUNCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	