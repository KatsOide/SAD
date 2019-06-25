from Tkinter import *

class _setit:
	def __init__(self, var, value):
		self.__value = value
		self.__var = var
	
	def __call__(self, *args):
		self.__var.set(self.__value)

class Optionmenu(Menubutton):
	def __init__(self, master, cnf={}, **kw):
		Widget.__init__(self, master, "menubutton", cnf, kw)
		self.widgetName = 'tk_optionMenu'
	
	def set_variable(self, variable=None, *values):
		kw = {"borderwidth": 2, "textvariable": variable,
		      "indicatoron": 1, "relief": RAISED, "anchor": "c",
		      "highlightthickness": 2}
#		Widget.__init__(self, master, "menubutton", kw)
		self.config(kw)
		menu = self.__menu = Menu(self, name="menu", tearoff=0)
		self.menuname = menu._w
		for v in values:
			menu.add_command(label=v, command=_setit(variable, v))
		self["menu"] = menu
	
	def __getitem__(self, name):
		if name == 'menu':
			return self.__menu
		return Widget.__getitem__(self, name)
	
	def destroy(self):
		Menubutton.destroy(self)
		self.__menu = None
	
	def set_menu(self, menu):
		self.menuname = menu._w
