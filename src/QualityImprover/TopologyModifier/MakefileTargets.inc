all_headers: $(TOPOMODIFHEADERS)
all_objects: $(TOPOMODIFOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(TOPOMODIFDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(TOPOMODIFDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(TOPOMODIFDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	