#!python
#-*- coding:utf-8 -*-
"""
converting a text surrounded by EquAlign[Text[]] in Packages/HelpMessage.n 
using pandoc program.
current version of pandoc converts "\begin{align}" in an input text to
"\[\begin{aligned". This causes NO numbered equations. sed is used to
fix this problem.
"""
from bs4 import BeautifulSoup as BS, NavigableString as NS, Tag
import urllib, urllib2
import datetime,time
import json
import os


soup=BS(open("../script/SADHelp-MJ.html"),'lxml')
# when this program is called in help2HTP-MV.sad, "cwd" is "Documents" directory.

t2h=soup.find_all('div',{'class':'Text2HTML'})

for i,e in enumerate(t2h):
    # pandoc replace \begin{align} to \[\begin{aligned}. sed is used to recover this behaviour.
    pin,pout=os.popen2(u"/usr/local/bin/pandoc --mathjax --from=latex --to=html|sed s/aligned\}/align\}/g ")
    pin.write(e.string)
    pin.close()
    o=pout.read()
    pout.close()
    b=BS(o,'lxml')
    bl=b.body.contents
    #print i,j,len(b.body.contents),type(c),map(type, b.body.contents)
    #print i,j,len(b.body.contents),b.body.contents
    e.contents[0].replace_with(bl[0]) # the elements in b.body.contents is cosumed by this process
    #print i,"len:",len(bl), type(e),"E:", e
    #print "Contents:", bl
    #print "BL:",bl
    k=1
    while b.body.contents:
        bc=b.body.contents[0]
        #print i, k, type(e), len(e.contents), type(bc), bc
        e.insert(k, bc)
        k +=1
    #print i, "E: ", e
    #print i,"t2h: ", t2h[i]
    
f=open("../Documents/SADHelp-MJ.html","w")
f.write(str(soup))
f.close()

