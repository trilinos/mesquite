#include "FileTokenizer.hpp"

#include <string>
#include <string.h>
using std::string;

namespace Mesquite
{

FileTokenizer::FileTokenizer( FILE* file_ptr )
  : filePtr( file_ptr ),
    nextToken( buffer ),
    bufferEnd( buffer ),
    lineNumber( 1 ),
    lastChar( '\0' )
  {}
  
FileTokenizer::~FileTokenizer() 
  { fclose( filePtr ); }

bool FileTokenizer::eof() const
  { return nextToken == bufferEnd && feof(filePtr); }

const char* FileTokenizer::get_string( MsqError& err )
{
    // If the whitepsace character marking the end of the
    // last token was a newline, increment the line count.
  if (lastChar == '\n')
    ++lineNumber;
  
    // Loop until either found the start of a token to return or have
    // reached the end of the file.
  for (;;)
  {
      // If the buffer is empty, read more.
    if (nextToken == bufferEnd)
    {
      size_t count = fread( buffer, 1, sizeof(buffer) - 1, filePtr );
      if (!count)
      {
        if (feof(filePtr))
          err.set_msg( "File truncated.\n");
        else
          err.set_msg( "I/O error.\n");
        return NULL;
      }
      
      nextToken = buffer;
      bufferEnd = buffer + count;
    }
    
      // If the current character is not a space, we've found a token.
    if (!isspace(*nextToken))
      break;
      
      // If the current space character is a newline,
      // increment the line number count.
    if (*nextToken == '\n')
      ++lineNumber;
    ++nextToken;
  }
  
    // Store the start of the token in "result" and
    // advance "nextToken" to one past the end of the
    // token.
  char* result = nextToken;
  while (nextToken != bufferEnd && !isspace(*nextToken))
    ++nextToken;
  
    // If we have reached the end of the buffer without finding
    // a whitespace character terminating the token, we need to
    // read more from the file.  Only try once.  If the token is
    // too large to fit in the buffer, give up.
  if (nextToken == bufferEnd)
  {
      // Shift the (possibly) partial token to the start of the buffer.
    size_t remaining = bufferEnd - result;
    memmove( buffer, result, remaining );
    result = buffer;
    nextToken = result + remaining;
    
      // Fill the remainder of the buffer after the token.
    size_t count = fread( nextToken, 1, sizeof(buffer) - remaining - 1, filePtr );
    if (!count && !feof(filePtr))
    {
      err.set_msg( "I/O error.\n");
      return NULL;
    }
    bufferEnd = nextToken + count;
    
      // Continue to advance nextToken until we find the space
      // terminating the token.
    while (nextToken != bufferEnd && !isspace(*nextToken))
      ++nextToken;
  
    if (nextToken == bufferEnd) // EOF
    {
      *bufferEnd = '\0';
      ++bufferEnd;
    }
  }
  
    // Save terminating whitespace character (or NULL char if EOF).
  lastChar = *nextToken;
    // Put null in buffer to mark end of current token.
  *nextToken = '\0';
    // Advance nextToken to the next character to search next time.
  ++nextToken;
  return result;
}

bool FileTokenizer::get_double_internal( double& result, MsqError& err )
{
    // Get a token
  const char *token_end, *token = get_string( err );
  if (!token)
    return false;
  
    // Parse token as double
  result = strtod( token, (char**)&token_end );

    // If the one past the last char read by strtod is
    // not the NULL character terminating the string,
    // then parse failed.
  if (*token_end)
  {
    char buffer[16];
    sprintf( buffer, "%d", line_number() );
    string message( "Syntax error at line " );
    message += buffer;
    message += ": expected number, got \"";
    message += token;
    message += "\"";
    err.set_msg( message );
    return false;
  }
  
  return true;
}

bool FileTokenizer::get_float_internal( float& result, MsqError& err )
{
  double d;
  if (!get_double_internal( d, err ))
    return false;
  
  result = (float)d;
  if (d != (double)result)
  {
    set_error( err, "Numeric overflow" );
    return false;
  }
  
  return true;
}

bool FileTokenizer::get_long_int_internal( long& result, MsqError& err )
{
    // Get a token
  const char *token_end, *token = get_string( err );
  if (!token)
    return false;
  
    // Parse token as long
  result = strtol( token, (char**)&token_end, 0 );

    // If the one past the last char read by strtol is
    // not the NULL character terminating the string,
    // then parse failed.
  if (*token_end)
  {
    char buffer[16];
    sprintf( buffer, "%d", line_number() );
    string message( "Syntax error at line " );
    message += buffer;
    message += ": expected integer number, got \"";
    message += token;
    message += "\"";
    err.set_msg( message );
    return false;
  }

  return true;
}

bool FileTokenizer::get_byte_internal( unsigned char& result, MsqError& err )
{
  long i;
  if (!get_long_int_internal( i, err ))
    return false;
  
  result = (unsigned char)i;
  if (i != (long)result)
  {
    set_error( err, "Numeric overflow" );
    return false;
  }
  
  return true;
}

bool FileTokenizer::get_short_int_internal( short& result, MsqError& err )
{
  long i;
  if (!get_long_int_internal( i, err ))
    return false;
  
  result = (short)i;
  if (i != (long)result)
  {
    set_error( err, "Numeric overflow" );
    return false;
  }
  
  return true;
}

bool FileTokenizer::get_integer_internal( int& result, MsqError& err )
{
  long i;
  if (!get_long_int_internal( i, err ))
    return false;
  
  result = (int)i;
  if (i != (long)result)
  {
    set_error( err, "Numeric overflow" );
    return false;
  }
  
  return true;
}

bool FileTokenizer::get_boolean_internal( bool& result, MsqError& err )
{
    // Get a token
  const char *token = get_string( err );
  if (!token)
    return false;
  
  if (token[1] || (token[0] != '0' && token[1] != '1'))
  {
    char buffer[16];
    sprintf( buffer, "%d", line_number() );
    string message( "Syntax error at line " );
    message += buffer;
    message += ": expected integer number, got \"";
    message += token;
    message += "\"";
    err.set_msg( message );
    return false;
  }

  result = token[0] == '1';
  return true;
}

bool FileTokenizer::get_floats( size_t count, float* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_float_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

bool FileTokenizer::get_doubles( size_t count, double* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_double_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

bool FileTokenizer::get_bytes( size_t count, unsigned char* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_byte_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

bool FileTokenizer::get_short_ints( size_t count, short* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_short_int_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}


bool FileTokenizer::get_integers( size_t count, int* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_integer_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

bool FileTokenizer::get_long_ints( size_t count, long* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_long_int_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

bool FileTokenizer::get_booleans( size_t count, bool* array, MsqError& err )
{
  for (size_t i = 0; i < count; ++i)
  {
    if (!get_boolean_internal( *array, err ))
      return false;
    ++array;
  }
  return true;
}

void FileTokenizer::unget_token()
{
  if (nextToken - buffer < 2)
    return;
  
  --nextToken;
  *nextToken = lastChar;
  --nextToken;
  while (nextToken > buffer && *nextToken)
    --nextToken;
    
  if (!*nextToken)
    ++nextToken;
    
  lastChar = '\0';
}

bool FileTokenizer::match_token( const char* str, MsqError& err )
{
    // Get a token
  const char *token = get_string( err );
  if (!token)
    return false;

    // Check if it matches
  if (0 == strcmp( token, str ))
    return true;
  
    // Construct error message
  string message( "Parsing error at line " );
  char lineno[16];
  sprintf( lineno, "%d", line_number() );
  message += lineno;
  message += ": expected \"";
  message += str;
  message += "\", got \"";
  message += token;
  message += "\"";
  err.set_msg( message );
  return false;
}  // namespace Mesquite


int FileTokenizer::match_token( const char* const* list, MsqError& err )
{
    // Get a token
  const char *token = get_string( err );
  if (!token)
    return false;

    // Check if it matches any input string
  const char* const* ptr;
  for (ptr = list; *ptr; ++ptr)
    if (0 == strcmp( token, *ptr ))
      return ptr - list + 1;
  
    // No match, constuct error message
  string message( "Parsing error at line " );
  char lineno[16];
  sprintf( lineno, "%d", line_number() );
  message += lineno;
  message += ": expected one of {";
  for (ptr = list; *ptr; ++ptr)
  {
    message += " ";
    message += *ptr;
  }
  message += " } got \"";
  message += token;
  message += "\"";
  err.set_msg( message );
  return false;
}

bool FileTokenizer::get_newline( MsqError& err )
{
  if (lastChar == '\n')
  {
    lastChar = ' ';
    ++lineNumber;
    return true;
  }
  
    // Loop until either we a) find a newline, b) find a non-whitespace
    // character or c) reach the end of the file.
  for (;;)
  {
      // If the buffer is empty, read more.
    if (nextToken == bufferEnd)
    {
      size_t count = fread( buffer, 1, sizeof(buffer), filePtr );
      if (!count)
      {
        if (eof())
          err.set_msg( "File truncated.\n");
        else
          err.set_msg( "I/O error.\n");
        return false;
      }
      
      nextToken = buffer;
      bufferEnd = buffer + count;
    }
    
      // If the current character is not a space, the we've failed.
    if (!isspace(*nextToken))
    {
      set_error( err, "Expected newline" );
      return false;
    }
      
      // If the current space character is a newline,
      // increment the line number count.
    if (*nextToken == '\n')
    {
      ++lineNumber;
      ++nextToken;
      lastChar = ' ';
      return true;
    }
    ++nextToken;
  }
  
    // should never reach this
  return false;
}


void FileTokenizer::set_error( MsqError& err, const char* prefix )
{
  char line_str[16];
  sprintf( line_str, "%d", line_number() );
  std::string msg( prefix );
  msg += " at line ";
  msg += line_str;
  err.set_msg( msg );
}

}  // namespace Mesquite

