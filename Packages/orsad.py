from kekbdatabase import *

def pr(x):print x

def printor(x):
        a=kekbdb(x)
        print("****description****")
        map(pr,a.description)
        print("****data****")
        map(pr,a)
        print("****end****")

