all_headers: $(QMUNTANGLEHEADERS)
all_objects: $(QMUNTANGLEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMUNTANGLEDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QMUNTANGLEDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QMUNTANGLEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	