all_headers: $(QUALMETRICHEADERS)
all_objects: $(QUALMETRICOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALMETRICDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QUALMETRICDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QUALMETRICDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	