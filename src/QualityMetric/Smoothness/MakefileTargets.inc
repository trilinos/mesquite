all_headers: $(QMSMOOTHNESSHEADERS)
all_objects: $(QMSMOOTHNESSOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMSMOOTHNESSDIR)/%.hpp
	ln -sf ../$< $@

$(localincludedir)/%.hh : $(QMSMOOTHNESSDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QMSMOOTHNESSDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	