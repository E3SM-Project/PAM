
#ifndef _PARSER_H_
#define _PARSER_H_

#include "const.h"
#include <fstream>
#include <string>
#include "string.h"
#include <sstream>
#include <algorithm>
#include "Domain.h"
#include "FileIO.h"
#include "Parallel.h"


void readParamsFile(std::string fNameIn, Domain &dom, Parallel &par, FileIO &io);

void handleTimeMethod(std::string &str, std::string &fNameIn);

void handleEqnSet(std::string &str, std::string &fNameIn, Domain &dom);

void handleDataInit(std::string &str, std::string &fNameIn, Domain &dom);


#endif
