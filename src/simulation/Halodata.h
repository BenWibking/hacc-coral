//-*-C++-*-/////////////////////////////////////////////////////////////
// $Id: Halodata.h,v 1.1 2010/03/11 04:14:06 pope Exp $
////////////////////////////////////////////////////////////////////////

#ifndef HALODATA_H
#define HALODATA_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdexcept>



class Halodata {

public:
  
  explicit Halodata( const char *file ) :
    m_bb(-1.0),
    m_minmass(-1.0),
    m_nmin(-1),
    m_oL(-1.0),
    m_Nskewers(-1),
    m_Npixels(-1),
    m_h(-1.0),
    m_Nmesh(-1)
  {
    std::ifstream myfile;

    myfile.open(file);
    if (! myfile.is_open()) {
      std::ostringstream ost;
      ost << "Halodata: cannot open '"
	  << file
	  << "'";
      throw std::runtime_error( ost.str() ); 
    }

    read_item( myfile, &m_bb );
    read_item( myfile, &m_minmass );
    read_item( myfile, &m_nmin );
    read_item( myfile, &m_oL );
    read_item( myfile, &m_Nskewers );
    read_item( myfile, &m_Npixels );
    read_item( myfile, &m_h );
    read_item( myfile, &m_Nmesh );

  }

  float bb() const { return m_bb; }
  float minmass() const { return m_minmass; }
  int nmin() const { return m_nmin; }
  float oL() const { return m_oL; }

  int Nskewers() const { return m_Nskewers; }
  int Npixels() const { return m_Npixels; }
  float h() const { return m_h; }
  int Nmesh() const { return m_Nmesh; }

private:

  Halodata();
  Halodata( const Halodata& );
  Halodata& operator = ( const Halodata& );

  template< typename U > void 
  read_item( std::ifstream& infile, U* item )
  {
    infile >> *item;
    std::string tmp_str;
    getline( infile, tmp_str );
    //if ((! infile) && (! infile.eof())) {
    if ( infile.eof() ) {
      throw std::runtime_error( "Halodata: error reading input file" );
    }
  }

  //halo parameters
  float m_bb;
  float m_minmass;
  int m_nmin;
  float m_oL;

  //skewer parameters
  int m_Nskewers;
  int m_Npixels;
  float m_h;
  int m_Nmesh;
};



#endif
