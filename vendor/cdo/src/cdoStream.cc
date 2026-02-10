/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdoStream.h"

static int nextCdoStreamID = 1;

CdoStream::~CdoStream() {}
CdoStream::CdoStream() : m_cdoStreamID(nextCdoStreamID++) {}

int
CdoStream::getTsID()
{
  return m_tsID;
}
