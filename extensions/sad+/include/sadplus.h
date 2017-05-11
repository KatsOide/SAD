#ifndef WALNUT_H
#define WALNUT_H

#include <cstring>
using std::strlen;
using std::strstr;

#include <dacpp.h>
#include <c_da.h>
#include <element.h>
#include <map_double.h>
#include <map_da.h>
#include <map_c_da.h>
#include <matrix.h>
#include <lin_map.h>
#include <lie_da.h>

void resource_files(int,char**,char*);
void user_interface(Accelerator);
int command_intp(char*,char*,char**); 
int get_com_id(char*);
void rm_space(char*);
void Tracking(Accelerator);
void Tracking(map_da);
void DAmapping(Accelerator);
void hamilton(const map_da&);

//int Stack_Search(char*,Stack_damap&,Stack_da&a,Stack_map_double&);
//void Stack_Remove(char*,int,Stack_damap&,Stack_da&,Stack_map_double&);
//void set_return(char*,Stack_damap&,Stack_da&,Stack_map_double&);
enum { TRACK=0, DAMAP,   HELP,   EXIT,  SYMPLECTIC,
       FAC1,    FAC2,    FAC3,   NFORM,    DAINIT,
       PRINT,   WRITE,   READ,   AGAIN,    VARTYPE,
        EQUAL ,    HAMILTON};
const char *command_list[]=
    {"Track",  "DAmap", "Help",  "Exit",  "Is_symplectic",
     "Fac_1",  "Fac_2", "Fac_3", "NormalForm","DAinit",
     "Print",  "Write",  "Read",  "Again",   "VarType",
     "Equal", "Hamilton",     ""};
enum { NOTHING=0, IN_DAMAP, IN_DA, IN_TRACK, IN_DOUBLE };


inline char* com_parm(char* s,int com_id) 
{
   return strstr(s,command_list[com_id])+strlen(command_list[com_id])+1;
}

#endif
