"""sad.py: calling sad functions from Python interpreter."""

import _sad

eval=_sad.eval
evalcommand=_sad.evalcommand
interp=_sad.interp

# please add new functionarities below. NY

def event(event,com,tag): _sad.eval(
  "Block[{$Event={Widget:>"+event.widget.widgetName+
  ",Tag->"+retv(tag)+
  ",Type->EventType["+`event.type`+
  "],X->"+`event.x`+
  ",Y->"+`event.y`+
  ",XRoot->"+`event.x_root`+
  ",YRoot->"+`event.y_root`+
  ",Height->"+`event.height`+
  ",Width->"+`event.width`+
  ",Char->"+`event.char`+
  ",KeySym->"+`event.keysym`+
#  ",Focus->"+`event.focus`+
#  ",SendEvent->"+`event.send_event`+
  ",KeyCode->"+`event.keycode`+
  ",State->"+`event.state`+
  ",KeySymNum->"+`event.keysym_num`+
  ",Time->"+`event.time`+
  "}},If[Py$EchoValue,Print['$Event=',$Event]];"+com+"]")

def ret(var):
  _sad.eval("Py$Return="+retv(var))

def retv(var):
  if type(var) == type((1,2)):
    return "{"+`var`[1:-1]+"}"
  elif type(var) == type('a'):
    return "\""+var+"\""
  elif type(var) == type([1,2]):
    return "{"+`var`[1:-1]+"}"
  else:
    return `var`
