all_headers: $(STEEPDESCHEADERS)
all_objects: $(STEEPDESCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(STEEPDESCDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(STEEPDESCDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(STEEPDESCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	