all_headers: $(CONTROLHEADERS)
all_objects: $(CONTROLOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(CONTROLDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(CONTROLDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(CONTROLDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	