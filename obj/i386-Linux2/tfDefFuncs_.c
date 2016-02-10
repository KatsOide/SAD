/* tfDefFunc_.c generated at 2010-12-28 17:08:30 -0800 */
extern int init_framework_Feature(void);
extern int init_framework_Random(void);
extern int init_framework_RandomPlugin_SAD(void);
extern int init_framework_FuncTBL(void);
extern int init_framework_DynlMI(void);
extern int sadDefFunc_Process(void);
extern int sadDefFunc_Crypt(void);
extern int sadDefFunc_FileIO(void);
extern int sadDefFunc_FeatureQ(void);
extern int sadDefFunc_DateTAI(void);
extern int sadDefFunc_Random(void);
extern int sadDefFunc_DynamicLink(void);
extern int sadDefFunc_NetworkIO(void);
extern int sadDefFunc_NetSemaphore(void);
extern int sadDefFunc_RenumberElement(void);
extern int sadDefFunc_Xlib(void);
extern int sadDefFunc_TkInter(void);
extern int sadDefFunc_BuildInfo(void);

int tfdeffuncs_(void) {
/* Call init_frameworks */
  init_framework_Feature();
  init_framework_Random();
  init_framework_RandomPlugin_SAD();
  init_framework_FuncTBL();
  init_framework_DynlMI();

/* Call sadDefFuncs */
  sadDefFunc_Process();
  sadDefFunc_Crypt();
  sadDefFunc_FileIO();
  sadDefFunc_FeatureQ();
  sadDefFunc_DateTAI();
  sadDefFunc_Random();
  sadDefFunc_DynamicLink();
  sadDefFunc_NetworkIO();
  sadDefFunc_NetSemaphore();
  sadDefFunc_RenumberElement();
  sadDefFunc_Xlib();
  sadDefFunc_TkInter();
  sadDefFunc_BuildInfo();

  return 0;
}

/* End of File */
