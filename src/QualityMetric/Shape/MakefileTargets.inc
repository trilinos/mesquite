all_headers: $(QMSHAPEHEADERS)
all_objects: $(QMSHAPEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMSHAPEDIR)/%.hpp
	ln -sf ../$< $@
	
$(localincludedir)/%.hh : $(QMSHAPEDIR)/%.hh
	ln -sf ../$< $@

$(localobjdir)/%.o : $(QMSHAPEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} ${CXXFLAGS} -c -o $@ $<

	