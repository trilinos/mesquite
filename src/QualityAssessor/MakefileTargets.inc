all_headers: $(QUALASSESSHEADERS)
all_objects: $(QUALASSESSOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALASSESSDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QUALASSESSDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QUALASSESSDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	