#include "MC3Options.h"

#include <cstdio>
#include <cstring>
using namespace std;

MC3Options::MC3Options(int argc, char *argv[]) :
  m_gridQ(0),
  m_restartDumpQ(0),
  m_pkDumpQ(0),
  m_aliveDumpQ(0),
  m_refreshQ(0),
  m_analysisQ(0),
  m_staticQ(0),
  m_lightconeQ(0),
  m_skewerQ(0),
  m_haloQ(0),
  m_initialAlltoallQ(false),
  m_interpQ(0),
  m_nInterp(0),
  m_polyQ(0),
  m_bladeQ(0),
  m_cmQ(0),
  m_everyCMQ(0),
  m_skipStreamQ(0),
  m_skipKickLRQ(0),
  m_skipKickSRQ(0),
  m_dontDropMemory(0),
  m_dontUseBigchunk(0),
  m_useFastTreeEval(0),
  m_useRCBTree(0),
  m_rcbTreeExtraLevels(0),
  m_rcbTreePPN(0),
  m_rcbTreeTaskPartMin(0),
  m_whiteNoiseInit(0),
  m_mpiio(0),
  m_topologyQ(0)
{
  m_exeName = argv[0];

  int gor;
  opterr = 0;
  while( (gor = getopt(argc, argv, "gMBmhzr:p:a:o:s:l:f:i:Pbce123FSRL:N:OT:wIt:") ) != -1)
    switch (gor) {
      case 'M':
        m_dontDropMemory = 1;
        break;
      case 'B':
        m_dontUseBigchunk = 1;
        break;
      case 'g':
	m_gridQ = 1;
	break;
      case 'm':
	m_initialAlltoallQ = true;
	break;
      case 'h':
	m_haloQ = 1;
	break;
      case 'z':
	m_skewerQ = 1;
	break;
      case 'r':
	m_restartDumpQ = 1;
	m_restartDumpName = optarg;
	break;
      case 'p':
	m_pkDumpQ = 1;
	m_pkDumpName = optarg;
	break;
      case 'a':
	m_aliveDumpQ = 1;
	m_aliveDumpName = optarg;
	break;
      case 'o':
	m_analysisQ = 1;
	m_analysisdatName = optarg;
	break;
      case 's':
	m_staticQ = 1;
	m_staticDumpName = optarg;
	break;
      case 'l':
	m_lightconeQ = 1;
	m_LCUpdateName = optarg;
	break;
      case 'f':
	m_refreshQ = 1;
	m_refreshName = optarg;
	break;
      case 'i':
	m_interpQ = 1;
	m_nInterp = atoi(optarg);
	break;
      case 'P':
	m_polyQ = 1;
	break;
      case 'b':
	m_bladeQ = 1;
	break;
      case 'c':
	m_cmQ = 1;
	break;
      case 'e':
	m_everyCMQ = 1;
	break;
      case '1':
	m_skipStreamQ = 1;
	break;
      case '2':
	m_skipKickLRQ = 1;
	break;
      case '3':
	m_skipKickSRQ = 1;
	break;
      case 'F':
	m_useFastTreeEval = 1;
	break;
      case 'S':
        if (m_useRCBTree == 0 || m_useRCBTree == 3) m_useRCBTree += 1;
        break;
      case 'R':
        if (m_useRCBTree == 0 || m_useRCBTree == 3) m_useRCBTree += 2; /* monopole */
        break;
      case 'L':
        m_rcbTreeExtraLevels = atoi(optarg);
        break;
      case 'N':
        m_rcbTreePPN = atoi(optarg);
        break;
      case 'T':
        m_rcbTreeTaskPartMin = atoi(optarg);
        break;
      case 'O':
        if (m_useRCBTree < 3) m_useRCBTree += 3;
        break;
      case 'w':
	m_whiteNoiseInit = 1;
	break;
      case 'I':
	m_mpiio = 1;
	break;
      case 't':
	m_topologyQ = 1;
	m_topologyString = optarg;
	break;
      case '?':
	if (isprint(optopt))
	  fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	else
	  fprintf(stderr,"Unknown option character `\\x%x'.\n", optopt);
	exit(-1);
      default:
	exit(-1);
    }
}
