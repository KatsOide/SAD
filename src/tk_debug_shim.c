// Tk 8.3 時代の内部デバッグ変数。現行 Tk には無いので自前で提供する。
// weak にしておくと、本物があればそちらが優先される（衝突回避）。
__attribute__((weak)) int bTkCanvTextDebug = 0;
