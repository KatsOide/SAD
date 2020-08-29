#include <sim/sad_api.h>
#include <sim/TFCODE.h>

#include <stdio.h>
#include <dbStaticLib.h>

static DBBASE *pdbbase = NULL;

void cdbstr_(int *nfunc, char *str, int *iret)
{
  switch (*nfunc) {
  case 1:
	/* dbReadDatabase (not used now) */
	*iret = dbReadDatabase(&pdbbase, str, NULL, NULL);
	return;
  case 2:
	/* dbWriteRecord */
	*iret = dbWriteRecord(pdbbase, str, NULL, 0);
	return;
  case 8:
	/* dbPath */
	*iret = dbPath(pdbbase, str);
	return;
  case 9:
	/* dbAddPath */
	*iret = dbAddPath(pdbbase, str);
	return;
  case 66:
	/* dbDumpRecord */
	dbDumpRecord(pdbbase, *str?str:NULL, 0);
	*iret = 0;
	return;
  }
}

void cdbstrx3_(int *nfunc, char *str1, char *str2, char *str3, int *iret)
{
  switch (*nfunc) {
  case 1:
	/* dbReadDatabase */
	*iret = dbReadDatabase(&pdbbase, str1, str2, *str3?str3:NULL);
	return;
  }
}

void cdb2int_(int *nfunc, int *nret, int *iret)
{
  DBENTRY *dbe;

  switch (*nfunc) {
  case 12:
	/* dbFirstRecordType */
	dbe = dbAllocEntry(pdbbase);
	if (dbe) {
	  *iret = dbFirstRecordType(dbe);
	  *nret = (int)dbe;
	} else
	  *nret = 0;
	return;
  }
}

void cdbint2str_(integer8 *kx, integer4 *nfunc, integer4 *iarg)
{
  char *str = NULL;

  switch(*nfunc) {
  case 14: /* dbGetRecordTypeName */
    str = dbGetRecordTypeName((DBENTRY *)*iarg);
    break;
  case 19: /* dbGetFieldName */
    str = dbGetFieldName((DBENTRY *)*iarg);
    break;
  case 20: /* dbGetDefault */
    str = dbGetDefault((DBENTRY *)*iarg);
    break;
  case 40: /* dbGetString */
    str = dbGetString((DBENTRY *)*iarg);
    break;
  }
  if(str != NULL) {
    *kx = ktfstring + ktsalocb(-1, str);
  } else {
    *kx = ktfoper + mtfnull;
  }
}

void cdbint2int_(int *nfunc, int *iarg, int *nret, int *iret)
{
  DBENTRY *dbe;

  switch (*nfunc) {
  case 10:
	/* dbGetNRecordTypes */
	*nret = dbGetNRecordTypes((DBENTRY *)*iarg);
	return;
  case 13:
	/* dbNextRecordType */
	dbe = dbCopyEntry((DBENTRY *)*iarg);
	if (dbe) {
	  *iret = dbNextRecordType(dbe);
	  *nret = (int)dbe;
	} else
	  *nret = 0;
	return;
  case 18:
	/* dbGetFieldType */
	*nret = dbGetFieldType((DBENTRY *)*iarg);
	return;
  case 45:
	/* dbGetNMenuChoices */
	*nret = dbGetNMenuChoices((DBENTRY *)*iarg);
	return;
  }
}

void cdbintint2int_(int *nfunc, int *iarg1, int *iarg2, int *nret, int *iret)
{
  DBENTRY *dbe;

  switch (*nfunc) {
  case 15:
	/* dbGetNFields */
	*nret = dbGetNFields((DBENTRY *)*iarg1, *iarg2);
	return;
  case 16:
	/* dbFirstField */
	dbe = dbCopyEntry((DBENTRY *)*iarg1);
	if (dbe) {
	  *iret = dbFirstField(dbe, *iarg2);
	  *nret = (int)dbe;
	} else
	  *nret = 0;
	return;
  case 17:
	/* dbNextField */
	dbe = dbCopyEntry((DBENTRY *)*iarg1);
	if (dbe) {
	  *iret = dbNextField(dbe, *iarg2);
	  *nret = (int)dbe;
	} else
	  *nret = 0;
	return;
  }
}

void cdbint2strs_(integer4 *nfunc, integer4 *iarg, integer8 *ps, integer4 *n)
{
  char **strs = NULL;
  int i;

  switch(*nfunc) {
  case 46: /* dbGetMenuChoices */
    *n = dbGetNMenuChoices((DBENTRY *)*iarg);
    if(*n > 0)
      strs = dbGetMenuChoices((DBENTRY *)*iarg);
    break;
  }
  if(*n > 0 && strs != NULL) {
    for(i = 0; i < *n; i++) {
      ps[i] = ktfstring + ktsalocb(-1, strs[i]);
    }
  } else {
    *n = -1;
  }
}

void cdbintstr2int_(int *nfunc, int *iarg, char *str, int *nret, int *iret)
{
  DBENTRY *dbe;

  switch (*nfunc) {
  case 30:
	/* dbCreateRecord */
	dbe = dbCopyEntry((DBENTRY *)*iarg);
	if (dbe) {
	  *iret = dbCreateRecord(dbe, str);
	  *nret = (int)dbe;
	} else
	  *nret = 0;
	return;
  case 38:
	/* dbFindField */
	/* printf("dbFindField: %d [%s]\n",*iarg,str); */
	dbe = dbCopyEntry((DBENTRY *)*iarg);
	if (dbe) {
	  *iret = dbFindField(dbe, str);
	  /*printf("dbFindField: %d\n",*iret);*/
	  if (!*iret) {
		*nret = (int)dbe;
		return;
	  }
	}
	dbFreeEntry(dbe);
	*nret = 0;
	return;
  }
}

void cdbintstr_(int *nfunc, int *iarg, char *str, int *iret)
{
  switch (*nfunc) {
  case 41:
	/* dbPutString */
	*iret = dbPutString((DBENTRY *)*iarg, str);
	return;
  }
}

void cdbintstr2str_(integer8 *kx, integer4 *nfunc, integer4 *iarg, char *str)
{
  char *rstr = NULL;

  switch(*nfunc) {
  case 42: /* dbVerify */
    /*    rstr = dbVerify((DBENTRY *)*iarg, str); */
    break;
  }
  if(rstr) {
    *kx = ktfstring + ktsalocb(-1, rstr);
  } else {
    *kx = ktfoper + mtfnull;
  }
}
