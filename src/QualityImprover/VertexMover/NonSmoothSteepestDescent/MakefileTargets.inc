all_headers: $(NONSMOOTHSTEEPDESCHEADERS)
all_objects: $(NONSMOOTHSTEEPDESCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(NONSMOOTHSTEEPDESCDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(NONSMOOTHSTEEPDESCDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(NONSMOOTHSTEEPDESCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	