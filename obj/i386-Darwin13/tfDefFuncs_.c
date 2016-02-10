/* tfDefFunc_.c generated at 2014-09-11 11:18:19 +0900 */
extern int init_framework_FuncTBL(void);
extern int init_framework_DynlMI(void);
extern int init_framework_Feature(void);
extern int init_framework_BuildInfo(void);
extern int init_framework_Random(void);
extern int init_framework_RandomPlugin_SAD(void);
extern int sadDefFunc_DynamicLink(void);
extern int sadDefFunc_FeatureQ(void);
extern int sadDefFunc_BuildInfo(void);
extern int sadDefFunc_Process(void);
extern int sadDefFunc_Crypt(void);
extern int sadDefFunc_FileIO(void);
extern int sadDefFunc_NetworkIO(void);
extern int sadDefFunc_NetSemaphore(void);
extern int sadDefFunc_TkInter(void);
extern int sadDefFunc_Xlib(void);
extern int sadDefFunc_DateTAI(void);
extern int sadDefFunc_Random(void);
extern int sadDefFunc_RenumberElement(void);

int tfdeffuncs_(void) {
/* Call init_frameworks */
  init_framework_FuncTBL();
  init_framework_DynlMI();
  init_framework_Feature();
  init_framework_BuildInfo();
  init_framework_Random();
  init_framework_RandomPlugin_SAD();

/* Call sadDefFuncs */
  sadDefFunc_DynamicLink();
  sadDefFunc_FeatureQ();
  sadDefFunc_BuildInfo();
  sadDefFunc_Process();
  sadDefFunc_Crypt();
  sadDefFunc_FileIO();
  sadDefFunc_NetworkIO();
  sadDefFunc_NetSemaphore();
  sadDefFunc_TkInter();
  sadDefFunc_Xlib();
  sadDefFunc_DateTAI();
  sadDefFunc_Random();
  sadDefFunc_RenumberElement();

  return 0;
}

/* End of File */
