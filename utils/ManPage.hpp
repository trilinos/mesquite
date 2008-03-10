/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ManPage.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MAN_PAGE_HPP
#define MSQ_MAN_PAGE_HPP

#include "Mesquite.hpp"
#include <iostream>
#include <string>

class ManPage
{
public:
  static msq_stdio::ostream& begin_bold( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".B" << msq_stdio::endl; }
  static msq_stdio::ostream& end_bold( msq_stdio::ostream& str )
    { return str << msq_stdio::endl; }
  static msq_stdio::ostream& bold( msq_stdio::ostream& str, const msq_std::string& s )
    { return end_bold( begin_bold(str) << s ); }

  static msq_stdio::ostream& begin_italic( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".I" << msq_stdio::endl; }
  static msq_stdio::ostream& end_italic( msq_stdio::ostream& str )
    { return str << msq_stdio::endl; }
  static msq_stdio::ostream& italic( msq_stdio::ostream& str, const msq_std::string& s )
    { return end_italic( begin_italic(str) << s ); }
    
  static msq_stdio::ostream& begin_section( msq_stdio::ostream& str, const msq_std::string& name )
    { return str << msq_stdio::endl << ".SH " << name << msq_stdio::endl; }
    
  static msq_stdio::ostream& begin_subsection( msq_stdio::ostream& str, const msq_std::string& name )
    { return str << msq_stdio::endl << ".SS " << name << msq_stdio::endl; }
    
  static msq_stdio::ostream& begin_paragraph( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".P " << msq_stdio::endl; }
    
  static msq_stdio::ostream& begin_hanging_paragraph( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".HP " << msq_stdio::endl; }
    
  static msq_stdio::ostream& begin_indent( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".RS " << msq_stdio::endl; }
  static msq_stdio::ostream& end_indent( msq_stdio::ostream& str )
    { return str << msq_stdio::endl << ".RE " << msq_stdio::endl; }
    
  static msq_stdio::ostream& begin_manpage( msq_stdio::ostream& str, const msq_std::string& name, int section )
    { return str << msq_stdio::endl << ".TH " << name << " " << section << msq_stdio::endl; }
    
  static msq_stdio::ostream& write_text( msq_stdio::ostream& str, bool hanging_indent, const msq_std::string& text );
};

#endif
